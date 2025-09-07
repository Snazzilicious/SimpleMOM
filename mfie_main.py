
import numpy as np
import sys
sys.path.insert(0, "/home/ianh/Documents/irad_ianh/empy")
import MeshUtils
import EM_Utils
import GeneralUtils

GeneralUtils.stdout("Beginning MOM-MFIE")

cmdline = GeneralUtils.get_standard_arg_parser()
args = cmdline.parse_args()


# load mesh
mesh_filename = args.mesh if args.mesh is not None else "/home/ianh/Documents/Aurora/Aurora-8.0.1-LINUX/tutorial/plane_wave/smooth_bumblebee.obj"
vertices, faces, groups, gnames = MeshUtils.load_obj( mesh_filename )
#vertices, faces = MeshUtils.load_stars( "/home/ianh/Documents/STARS/trunk/Validation/STARSruns/Pipe1/pipe1_full_m.facet.STARS.20200518194220" )
normals, v1s,v2s, areas = MeshUtils.local_basis_areas( vertices, faces )
centroids = MeshUtils.get_centroids( vertices, faces )
n_faces = len(faces)
GeneralUtils.stdout("Mesh Loaded")


# Excitation params
freq = args.frequency if args.frequency is not None else 1e9
w = EM_Utils.omega_from_Hz( freq )
k = EM_Utils.k_from_omega( w )

nrhs=1
z_ang = np.radians( 90.0 )
x_ang = np.radians( 0.0 )
p = -EM_Utils.x_angle_z_angle_to_dir( x_ang, z_ang )
H = -EM_Utils.x_angle_z_angle_to_dir( x_ang+np.pi, (np.pi/2) - z_ang )
incident = EM_Utils.planewave( k, p, H, centroids )



# construct impedance matrix
M = np.identity( 2*n_faces, dtype=np.complex128 )
for i in range(n_faces):
	g = EM_Utils.gradx_G( centroids[i], centroids, k )
	cross1 = np.cross( normals[i], np.cross( v1s, g ) )
	cross2 = np.cross( normals[i], np.cross( v2s, g ) )
	
	row11 = 2 * cross1.dot( v1s[i] ) * areas
	row12 = 2 * cross2.dot( v1s[i] ) * areas
	row21 = 2 * cross1.dot( v2s[i] ) * areas
	row22 = 2 * cross2.dot( v2s[i] ) * areas
	
	M[i,:n_faces] -= row11
	M[i,n_faces:] -= row12
	M[i+n_faces,:n_faces] -= row21
	M[i+n_faces,n_faces:] -= row22

# construct RHS
b = np.zeros( (2*n_faces,nrhs), dtype=np.complex128 )
for j in range(nrhs):
	n_cross_H = np.cross( normals, incident[:, 3*j:3*j+3] )
	
	b[:n_faces,j] = 2 * np.sum( v1s * n_cross_H, axis=1 )
	b[n_faces:,j] = 2 * np.sum( v2s * n_cross_H, axis=1 )

GeneralUtils.stdout("Impedance matrix and RHS constructed")

sol = np.linalg.solve( M, b )
J = ( v1s.T * sol[:n_faces,0] + v2s.T * sol[n_faces:,0] ).T

from pathlib import Path
casename = args.output if args.output is not None else "mfie_{0:s}_{1:.0f}Hz".format( Path(mesh_filename).stem, freq )
np.save(casename+"_current.npy", J)

GeneralUtils.stdout("Solution computed")

# compute farfield
obs_x_angs = np.linspace( np.radians(-180), np.radians(180), 361 )
obs_z_angs = np.radians(90) * np.ones( 361 )
obs_angs = np.column_stack( (obs_x_angs,obs_z_angs) )
ff = EM_Utils.bistatic_H_field( obs_angs, k, J, centroids, areas )

ff_x, ff_z = EM_Utils.projected_farfield( obs_angs, ff )
#plt.plot( np.degrees(obs_x_angs), EM_Utils.to_dBsm(ff_z) )

GeneralUtils.stdout("Farfield computed")


cell_results = { "H_inc_real" : np.real(incident), "H_inc_imag" : np.imag(incident) }
cell_results["J_real"] = np.real(J)
cell_results["J_imag"] = np.imag(J)
cell_results["areas"] = areas
cell_results["1/R"] = 1.0 / (np.linalg.norm( centroids[250]-centroids, axis=1 )+1e-15)

MeshUtils.write_vtk( casename+"_results.vtk", vertices, faces, cell_values=cell_results )











