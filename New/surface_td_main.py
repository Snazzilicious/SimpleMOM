
import numpy as np
import sys
sys.path.insert(0, "/home/ianh/Documents/miscstuff")
import Inputs
import Logging
import MeshUtils
import EM_Utils

Logging.stdout("Beginning Surface Time Domain")

cmdline = Inputs.get_standard_arg_parser()
args = cmdline.parse_args()


# load mesh
mesh_filename = args.mesh if args.mesh is not None else "/home/ianh/Documents/Aurora/Aurora-8.0.1-LINUX/tutorial/plane_wave/smooth_bumblebee.obj"
vertices, faces, groups, gnames = MeshUtils.load_obj( mesh_filename )
#vertices, faces = MeshUtils.load_stars( "/home/ianh/Documents/STARS/trunk/Validation/STARSruns/Pipe1/pipe1_full_m.facet.STARS.20200518194220" )
normals, v1s,v2s, areas = MeshUtils.local_basis_areas( vertices, faces )
centroids = MeshUtils.get_centroids( vertices, faces )
R = np.linalg.norm( centroids.reshape(centroids.shape+(1,)) - centroids.T.reshape((1,)+centroids.T.shape), axis=1 )
n_faces = len(faces)
Logging.stdout("Mesh Loaded")



# set initial values
rho = -np.ones( n_faces )
J1 = np.zeros( n_faces )
J2 = np.zeros( n_faces )

class F:
	def __init__(self):
		self.history = None
	
	def __call__( self, t, y ):
		
		rho_t = -Div( J )
		
		# TODO ( need to include contribution from positive charges (+1)e_p )
		E = ...
		B = ...
		
		J_t = ( e / m_e ) * ( rho*E + np.cross( J, B ) ) # TODO - Div( J * J/rho )
		J1_t = np.sum( v1s * J_t, axis=1 )
		J2_t = np.sum( v2s * J_t, axis=1 )
		
		return rho_t, J1_t, J2_t

# TODO
# Use FD and interpolation for time derivatives of history data

