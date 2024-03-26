
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0,epsilon_0,speed_of_light
from scipy.integrate import quad_vec
from numba import njit
from MeshUtils import loadVTK

# TODO
# verify derivates numerically
# compute pts and wts
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


# Set up Excitations
w = 3e9
k = w / speed_of_light

class PlaneWave:
	def __init__(self, propDir, polVec):
		self.prop = propDir.copy()
		self.polV = polVec.copy()
	
	def excitation(self, xyz ):
		return self.polV * np.exp( -1j * k * xyz.dot( self.prop ) )

excitations = [ PlaneWave( np.array([np.cos(th),np.sin(th),0]), np.array([0,0,1]) ) for th in np.linspace(0,np.pi,5) ]
'''
def integrand_dydxG( x,y, test_i,basis_j ):
	dx = x-y
	R = np.linalg.norm( dx )
	I = bases[test_i] @ ( np.outer(dx,dx) * (R**2*k**2 - 3j*R*k - 3) + (R**2 + 1j*R**3*k)*np.eye(3) ) @ bases[basis_j]
	return ( I / (1j * w * epsilon_0) ) * ( np.exp( -1j * k * R ) / ( 4 * np.pi * R**5 + 1e-9 ) )

'''

tst_pts = (2.0*vertices[facets[:,0],:] + vertices[facets[:,1],:] + vertices[facets[:,2],:])/4.0
src_pts = (vertices[facets[:,0],:] + vertices[facets[:,1],:] + 2.0*vertices[facets[:,2],:])/4.0

dx = tst_pts.reshape((-1,3,1)) - src_pts.reshape((-1,3,1)).T
R = np.linalg.norm( dx, axis=1 )

A_G = 1j * w * mu_0 * np.exp( -1j * k * R ) / ( 4 * np.pi * R ) # Plain ole Green's function
A_G *= np.outer( areas, areas )

A = np.zeros( (2*nFacets,2*nFacets), dtype=np.complex128 )

# Fill RHS
b = np.zeros( (2*nFacets,len(excitations)), dtype=np.complex128 )
for j,exc in enumerate(excitations):
	E_inc = exc.excitation( tst_pts )
	b[:,j] = -np.sum( bases * np.row_stack((E_inc,E_inc)), axis=1 )


# Factor and solve
sols = np.linalg.solve( A,b )


# post process
for j,obs_dir in enumerate(observations):
	phase = np.exp( -1j * k * src_pts @ obs_dir )
	farfield[j] = 1j * w * mu_0 * ( ( np.concatenate((phase,phase)) * sols[:,j] ) @ ( bases - np.outer( bases @ obs_dir, obs_dir ) ) ) # Farfield Green's function


# Plot results
#ax = plt.figure().add_subplot(projection='3d')
#ax.quiver( testPts[:,0],testPts[:,1],testPts[:,2], np.real(b[:,0]),np.real(b[:,1]),np.real(b[:,2]), length=0.1, normalize=True)
#plt.show()

plt.plot(np.linspace(0,2*np.pi,40),mags[0])
plt.show()

