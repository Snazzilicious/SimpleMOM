
import numpy as np
import sys
sys.path.insert(0, "/home/ianh/Documents/irad_ianh/empy")
import MeshUtils
import EM_Utils
import GeneralUtils
import POSBR

GeneralUtils.stdout("Beginning PO-SBR")

cmdline = GeneralUtils.get_standard_arg_parser()
cmdline.add_argument('--max_bounce', type=int, help='Maximum number of ray bounces')
cmdline.add_argument('--ray_spacing', type=float, help='Distance between rays')
args = cmdline.parse_args()


# load mesh
mesh_filename = args.mesh if args.mesh is not None else "/home/ian/Aurora-8.0.1-LINUX/tutorial/plane_wave/smooth_bumblebee.obj"
vertices, faces, groups, gnames = MeshUtils.load_obj( mesh_filename )
normals, v1s,v2s, areas = MeshUtils.local_basis_areas( vertices, faces )
centroids = MeshUtils.get_centroids( vertices, faces )
mesh_center = ( np.max(vertices,axis=0) + np.min(vertices,axis=0) )/2
mesh_diam = np.linalg.norm( np.max(vertices,axis=0) - np.min(vertices,axis=0) )
GeneralUtils.stdout("Mesh Loaded")


# Excitation params
freq = args.frequency if args.frequency is not None else 1e9
w = EM_Utils.omega_from_Hz( freq )
k = EM_Utils.k_from_omega( w )

z_ang = np.radians( 90.0 )
x_ang = np.radians( 0.0 )
p = -EM_Utils.x_angle_z_angle_to_dir( x_ang, z_ang )
H = -EM_Utils.x_angle_z_angle_to_dir( x_ang+np.pi, (np.pi/2) - z_ang )


# Create initial rays
ray_spacing = args.ray_spacing if args.ray_spacing is not None else 5e-2
rays_u = int( mesh_diam / ray_spacing ) +1
rays_v = int( mesh_diam / ray_spacing ) +1
n_init_rays = rays_u * rays_v
wave_front_center = mesh_center - mesh_diam * p
wave_front_u = np.array([ np.sin(x_ang), -np.cos(x_ang), 0.0 ])
wave_front_v = np.array([ -np.cos(x_ang)*np.cos(z_ang), -np.sin(x_ang)*np.cos(z_ang), np.sin(z_ang) ])
ray_dirs = np.tile( p, (n_init_rays,1) )
ray_orgs = np.tile( wave_front_center, (n_init_rays,1) )
for ind,org in enumerate(ray_orgs): # TODO vectorize
	i = ind // rays_u
	j = ind % rays_u
	u = ( i * ray_spacing ) - ( rays_u * ray_spacing )/2.0
	v = ( j * ray_spacing ) - ( rays_v * ray_spacing )/2.0
	ray_orgs[ind,:] += u*wave_front_u + v*wave_front_v
ray_pols = np.tile( H, (n_init_rays,1) )


GeneralUtils.stdout("Initial rays constructed")


# Shoot and bounce rays
n_rays = len(ray_dirs)
max_bounce_cnt = args.max_bounce if args.max_bounce is not None else 24
bounce_cnt = 0
tracer = POSBR.RayTracer( vertices, faces )
currents=[]
current_pts=[]
current_wts=[]
while bounce_cnt < max_bounce_cnt and n_rays > 0:
	# trace rays in queue to get collisions
	struck_face_inds = tracer.intersects_first( ray_orgs, ray_dirs )
	
	# preprocess hits 
	n_strikes = np.sum( struck_face_inds >= 0 )
	p = np.argsort( struck_face_inds )[-n_strikes:]
	struck_face_inds = struck_face_inds[p]
	ray_orgs = ray_orgs[p,:]
	ray_dirs = ray_dirs[p,:]
	ray_pols = ray_pols[p,:]
	
	# gather normals
	struck_normals = normals[struck_face_inds,:]
	dir_dot_normal = np.sum( ray_dirs*struck_normals, axis=1 )
	
	# compute strike locations
	t_fars = np.sum( (centroids[struck_face_inds] - ray_orgs) * struck_normals, axis=1 ) / dir_dot_normal
	ray_orgs = np.transpose( ray_orgs.T + t_fars * ray_dirs.T )
	
	# compute reflected direction
	ray_dirs = np.transpose( ray_dirs.T - 2 * dir_dot_normal * struck_normals.T )
	
	# compute polarization at strike point
	ray_pols = np.transpose( np.exp( -1j * k * t_fars ) * ray_pols.T )
	
	# compute currents
	J = 2 * np.cross( struck_normals, ray_pols )
	currents.append( J.copy() )
	current_pts.append( ray_orgs.copy() )
	current_wts.append( ray_spacing**2 / np.abs( dir_dot_normal ) )
	
	# compute reflected polarization
	ray_pols = np.transpose( np.cross( J, ray_dirs ).T / ( 2 * np.abs( dir_dot_normal ) ) )
	
	ray_orgs += 1e-3 * ray_dirs # XXX Prevents tracer from re-intersecting with origin face
	
	bounce_cnt += 1
	n_rays = len(ray_dirs)


J_final = np.vstack(currents)
J_final_pts = np.vstack(current_pts)
J_final_wts = np.concatenate(current_wts)

GeneralUtils.stdout("Solution computed")


# compute farfield
obs_x_angs = np.linspace( np.radians(-180), np.radians(180), 361 )
obs_z_angs = np.radians(90) * np.ones( 361 )
obs_angs = np.column_stack( (obs_x_angs,obs_z_angs) )
ff = EM_Utils.bistatic_H_field( obs_angs, k, J_final, J_final_pts, J_final_wts )

ff_x, ff_z = EM_Utils.projected_farfield( obs_angs, ff )
#plt.plot( np.degrees(obs_x_angs), EM_Utils.to_dBsm(ff_z) )

GeneralUtils.stdout("Farfield computed")
	
