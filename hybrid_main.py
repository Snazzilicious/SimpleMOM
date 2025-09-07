
import numpy as np
from scipy import sparse
from scipy.spatial import KDTree
import sys
sys.path.insert(0, "/home/ianh/Documents/irad_ianh/empy")
import MeshUtils
import EM_Utils
import GeneralUtils
import IPO

GeneralUtils.stdout("Beginning Hybrid IPO-MOM")

cmdline = GeneralUtils.get_standard_arg_parser()
args = cmdline.parse_args()


mesh_filename = args.mesh if args.mesh is not None else "/home/ianh/Documents/Aurora/Aurora-8.0.1-LINUX/tutorial/plane_wave/smooth_bumblebee.obj"
vertices, faces, groups, gnames = MeshUtils.load_obj( mesh_filename )
GeneralUtils.stdout("Mesh Loaded")
normals, v1s,v2s, areas = MeshUtils.local_basis_areas( vertices, faces )
centroids = MeshUtils.get_centroids( vertices, faces )
n_faces = len(faces)
GeneralUtils.stdout("Geometry calculations complete")


# Find MOM groups
nbr_map = MeshUtils.get_nbr_map( faces )
is_mom = MeshUtils.get_mom_faces( nbr_map, normals )


# Excitation params
freq = args.frequency if args.frequency is not None else 1e9
w = EM_Utils.omega_from_Hz( freq )
k = EM_Utils.k_from_omega( w )

# z and x angle positions of incident waves
incident_angs = [(90.0, 0.0)]
nrhs=len(incident_angs)
incident = np.zeros( (n_faces,3*nrhs), dtype=np.complex128 )
for j,(z_ang,x_ang) in enumerate(incident_angs):
	p = -EM_Utils.x_angle_z_angle_to_dir( x_ang, z_ang )
	H = -EM_Utils.x_angle_z_angle_to_dir( x_ang+np.pi, (np.pi/2) - z_ang )
	illuminated = IPO.is_illuminated_pw( p, normals )
	incident[:,3*j:3*j+3] = ( EM_Utils.planewave( k, p, H, centroids ).T * illuminated ).T


# construct impedance matrix
vis_map = IPO.get_recieve_map_optimized( centroids, normals )
# update visibility for MOM faces
kdt = KDTree(centroids)
nearest_neighbors = kdt.query_ball_point( centroids[is_mom], 2*(EM_Utils.c / freq) ) # includes self index TODO workers=-1
adtl_recvs = [[] for _ in vis_map]
for i,nbr_list in zip(np.where(is_mom)[0],nearest_neighbors):
	vis_map[i] = np.union1d( vis_map[i], nbr_list ) # update recieve list
	for j in nbr_list:
		adtl_recvs[j].append(i)
for i,recv_list in enumerate(adtl_recvs): # update radiate list
	if len(recv_list) > 0:
		vis_map[i] = np.union1d( vis_map[i], recv_list )


M = sparse.identity( 2*n_faces, dtype=np.complex128, format='lil' )
#M = sparse.lil_matrix( (2*n_faces,2*n_faces), dtype=np.complex128 )
for i in range(n_faces):
	if len(vis_map[i]) == 0:
		continue
	
	g = EM_Utils.gradx_G( centroids[i], centroids[vis_map[i]], k )
	cross1 = np.cross( normals[i], np.cross( v1s[vis_map[i]], g ) )
	cross2 = np.cross( normals[i], np.cross( v2s[vis_map[i]], g ) )
	
	A = 2 * areas[vis_map[i]]
	row11 = cross1.dot( v1s[i] ) * A
	row12 = cross2.dot( v1s[i] ) * A
	row21 = cross1.dot( v2s[i] ) * A
	row22 = cross2.dot( v2s[i] ) * A
	
	M[i,vis_map[i]] -= row11
	M[i,vis_map[i]+n_faces] -= row12
	M[i+n_faces,vis_map[i]] -= row21
	M[i+n_faces,vis_map[i]+n_faces] -= row22
print(M.nnz, "/", np.prod(M.shape))
M = M.tocsr()
	
# construct RHS / initial_current
b = np.zeros( (2*n_faces,nrhs), dtype=np.complex128 )
for j in range(nrhs):
	n_cross_H = np.cross( normals, incident[:, 3*j:3*j+3] )
	
	b[:n_faces,j] = 2 * np.sum( v1s * n_cross_H, axis=1 )
	b[n_faces:,j] = 2 * np.sum( v2s * n_cross_H, axis=1 )


GeneralUtils.stdout("Impedance matrix and RHS constructed")


sol = sparse.linalg.spsolve( M, b )
J = np.zeros( (n_faces,3*nrhs), dtype=np.complex128 )
for j in range(nrhs):
	J[:,3*j:3*j+3] = ( v1s.T * sol[:n_faces,j] + v2s.T * sol[n_faces:,j] ).T

from pathlib import Path
casename = args.output if args.output is not None else "ipomom_{0:s}_{1:.0f}Hz".format( Path(mesh_filename).stem, freq )
np.save(casename+"_current.npy", J)


GeneralUtils.stdout("Solution computed")


# compute farfield
obs_x_angs = np.linspace( np.radians(-180), np.radians(180), 361 )
obs_z_angs = np.radians(90) * np.ones( 361 )
obs_ang_sets = [np.column_stack(( obs_x_angs,obs_z_angs ))]
ffs = [{} for _ in range(nrhs)]
for j in range(nrhs):
	obs_angs = np.row_stack( obs_ang_sets[j] )
	ff = EM_Utils.bistatic_H_field( obs_angs, k, J[:,3*j:3*j+3], centroids, areas )
	ff_x, ff_z = EM_Utils.projected_farfield( obs_angs, ff )
	ffs[j]['x'] = ff_x
	ffs[j]['z'] = ff_z
	#plt.plot( np.degrees(obs_x_angs), EM_Utils.to_dBsm(ff_z) )


GeneralUtils.stdout("Farfield computed")


cell_results = { "H_inc_real" : np.real(incident), "H_inc_imag" : np.imag(incident) }
cell_results["J_real"] = np.real(J)
cell_results["J_imag"] = np.imag(J)
cell_results["areas"] = areas

MeshUtils.write_vtk( casename+"_result.vtk", vertices, faces, cell_values=cell_results )

""" Iterative solve
current_tol = 1e-2
max_its = 10

J0 = b

Js = [ J0 ]

err = np.inf
it=0
while err > current_tol and it < max_its :
	Js.append( b - M.dot( Js[-1] ) )
	err = np.linalg.norm( Js[-1] - Js[-2] ) / len(J0)
	it += 1

sol = sols[-1]
# TODO compute farfield
"""
