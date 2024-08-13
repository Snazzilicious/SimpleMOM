
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0,epsilon_0,speed_of_light
from MeshUtils import loadVTK

# TODO
# verify derivates numerically
# spot check matrix elements (and rhs)
# farfield
# find reference cases

# Load mesh
meshFileName = "SphereMesh.vtk"
vertices, facets = loadVTK( meshFileName )
nFacets = len(facets)
# Facet local vector basis
bases = np.zeros((2*nFacets,3))
bases[:nFacets,:] = vertices[facets[:,1],:] - vertices[facets[:,0],:]
bases[nFacets:,:] = vertices[facets[:,2],:] - vertices[facets[:,0],:]
areas = 0.5 * np.linalg.norm( np.cross( bases[:nFacets], bases[nFacets:] ), axis=1 )

# Set up Excitations
w = 3e9
k = w / speed_of_light

class PlaneWave:
	def __init__(self, propDir, polVec):
		self.prop = propDir.copy()
		self.polV = polVec.copy()
	
	def excitation(self, xyz ):
		return np.outer( np.exp( -1j * k * xyz.dot( self.prop ) ), self.polV )

excitations = [ PlaneWave( np.array([np.cos(th),np.sin(th),0]), np.array([0,0,1]) ) for th in np.linspace(0,np.pi,5) ]


# Fill Impedance Matrix
tst_pts = (2.0*vertices[facets[:,0],:] + vertices[facets[:,1],:] + vertices[facets[:,2],:])/4.0
src_pts = (vertices[facets[:,0],:] + vertices[facets[:,1],:] + 2.0*vertices[facets[:,2],:])/4.0

dx = tst_pts.reshape((-1,3,1)) - src_pts.reshape((-1,3,1)).T
R = np.linalg.norm( dx, axis=1 )

A_G = 1j * w * mu_0 * np.exp( -1j * k * R ) / ( 4 * np.pi * R ) # Plain ole Green's function

A_ddG_1 = R**2*k**2 - 3j*R*k - 3
A_ddG_2 = R**2 + 1j*R**3*k
A_ddG_3 = ( 1 / (1j * w * epsilon_0) ) * np.exp( -1j * k * R ) / ( 4 * np.pi * R**5 )

bTb = bases @ bases.T
bT_RRT_b_11 = np.einsum( "ijk,ij->ik", dx, bases[:nFacets] ) * np.einsum( "ijk,kj->ik", dx, bases[:nFacets] )
bT_RRT_b_12 = np.einsum( "ijk,ij->ik", dx, bases[:nFacets] ) * np.einsum( "ijk,kj->ik", dx, bases[nFacets:] )
bT_RRT_b_21 = np.einsum( "ijk,ij->ik", dx, bases[nFacets:] ) * np.einsum( "ijk,kj->ik", dx, bases[:nFacets] )
bT_RRT_b_22 = np.einsum( "ijk,ij->ik", dx, bases[nFacets:] ) * np.einsum( "ijk,kj->ik", dx, bases[nFacets:] )

weight = np.outer( areas, areas )

A = np.zeros( (2*nFacets,2*nFacets), dtype=np.complex128 )
A[:nFacets,:nFacets] = ( bTb[:nFacets,:nFacets] * ( A_G + A_ddG_2*A_ddG_3 ) + bT_RRT_b_11 * A_ddG_1*A_ddG_3 ) * weight
A[nFacets:,:nFacets] = ( bTb[nFacets:,:nFacets] * ( A_G + A_ddG_2*A_ddG_3 ) + bT_RRT_b_12 * A_ddG_1*A_ddG_3 ) * weight
A[:nFacets,nFacets:] = ( bTb[:nFacets,nFacets:] * ( A_G + A_ddG_2*A_ddG_3 ) + bT_RRT_b_21 * A_ddG_1*A_ddG_3 ) * weight
A[nFacets:,nFacets:] = ( bTb[nFacets:,nFacets:] * ( A_G + A_ddG_2*A_ddG_3 ) + bT_RRT_b_22 * A_ddG_1*A_ddG_3 ) * weight

# Fill RHS
b = np.zeros( (2*nFacets,len(excitations)), dtype=np.complex128 )
for j,exc in enumerate(excitations):
	E_inc = exc.excitation( tst_pts )
	b[:,j] = -np.sum( bases * np.row_stack((E_inc,E_inc)), axis=1 )


# Factor and solve
sols = np.linalg.solve( A,b )


# post process
for obs_dir,ex_ind in observations:
	phase = np.exp( -1j * k * src_pts @ obs_dir ) * areas
	obs.farfield = 1j * w * mu_0 * ( ( np.concatenate((phase,phase)) * sols[:,j] ) @ ( bases - np.outer( bases @ obs_dir, obs_dir ) ) ) # Farfield Green's function


plt.plot(np.linspace(0,2*np.pi,40),mags[0])
plt.show()

