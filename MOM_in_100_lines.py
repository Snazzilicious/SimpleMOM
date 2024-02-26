
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0,epsilon_0,speed_of_light
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
nVertices = len(vertices)
nFacets = len(facets)

# Construct vector basis functions
bases = np.zeros((2*nFacets,3))
bases[:nFacets] = vertices[facets[:,1],:] - vertices[facets[:,0],:]
bases[nFacets:] = vertices[facets[:,2],:] - vertices[facets[:,0],:]


# Set up Scenario(s)
w = 3e9
k = w / speed_of_light

class PlaneWave:
	def __init__(self, propDir, polVec, obs=None):
		self.prop = propDir.copy()
		self.polV = polVec.copy()
		
		if obs is None:
			self.observations = [-self.prop]
		else:
			self.observations = obs.copy()
	
	def excitation(self, xyz ):
		return self.polV * np.exp( 1j * k * self.prop.dot( xyz ) )

scenarios = [PlaneWave( np.array([1,0,0]), np.array([0,0,1]), [ np.array([np.cos(th),np.sin(th),0]) for th in np.linspace(0,np.pi,40) ] )]
scenarios.extend([ PlaneWave( np.array([np.cos(th),np.sin(th),0]), np.array([0,0,1]) ) for th in np.linspace(0,np.pi,5) ])


# Fill Matrix and RHS(s)
def G( x,y ):
	R = np.linalg.norm( x - y, axis=1 )
	return np.exp( 1j * k * R ) / ( 4 * np.pi * R + 1e-15 )

def gradxG( x,y ):
	dx = x-y
	R = np.linalg.norm( dx, axis=1 ).reshape(( x.shape[0], 1, y.shape[-1] ))
	return dx * ( 1j * k * R - 1 ) * np.exp( 1j * k * R ) / ( 4 * np.pi * R**3 + 1e-15 )

def gradygradxG( x,y ):
	dx = x-y
	R = np.linalg.norm( dx, axis=1 ).reshape(( x.shape[0], 1,1, y.shape[-1] ))
	prefix = np.einsum('ijk,ilk->ijlk', dx,dx) * (R**2*k**2 + 3j*R*k - 3) + (R**2 - 1j*R**3*k)*np.eye(3).reshape((x.shape[0],3,3,y.shape[-1]))
	return prefix * np.exp( 1j * k * R ) / ( 4 * np.pi * R**5 )


A = np.zeros((2*nFacets,2*nFacets),dtype=np.complex128)
b = np.zeros((2*nFacets,len(scenarios)),dtype=np.complex128)
for i in range(2*nFacets):
	for j in range(2*nFacets):
		facet_i = i % nFacets
		facet_j = j % nFacets
		# integral of  jwu G(x_i,x_j) vi^Tv_j + vi^T dxi dxj G(x_i,x_j) v_j / jwe
		I1 = 1j * w * mu_0 * (bases[i] @ bases[j]) * ( wts[facet_i] @ G(pts[facet_i],pts[facet_j]) @ wts[facet_j] )
		if facet_i == facet_j : # Singularity!
			# integral of -vi^T gradG vj^T nHat
			# TODO
			I2 = 
		else:
			I2 = bases[i] @ ( wts[facet_i] @ gradygradxG(pts[facet_i],pts[facet_j]) @ wts[facet_j] ) @ bases[j]
		I2 /= 1j * w * epsilon_0
		A[i,j] = I1 + I2
		
	for j in range(len(scenarios)):
		# integral of vi^T scenarios[j].excitation( xi )
		b[i,j] = wts[facet_i] @ scenarios[j].excitation(pts[facet_i]) @ bases[facet_i]

# TODO Replace slow double loop with these
A = 1j * w * mu_0 * bases @ bases.T
I_G = np.einsum('ijkl,ik,jl->ij', G(), wts,wts ) # TODO find most efficient operations & memory ordering
for i in range(2):
	for j in range(2):
		A[i*nFacets:(i+1)*nFacets,j*nFacets:(j+1)*nFacets] *= I_G

# Factor and solve
sols = np.linalg.solve( A,b )


# post process
mags = [[np.linalg.norm( np.exp( 1j * k * srcPts.dot( obs ) ) @ sols[:,3*j:3*(j+1)] ) for obs in scene.observations] for scene in scenarios]

# Plot results
#ax = plt.figure().add_subplot(projection='3d')
#ax.quiver( testPts[:,0],testPts[:,1],testPts[:,2], np.real(b[:,0]),np.real(b[:,1]),np.real(b[:,2]), length=0.1, normalize=True)
#plt.show()

plt.plot(np.linspace(0,2*np.pi,40),mags[0])
plt.show()

