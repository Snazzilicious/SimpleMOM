
import numpy as np
import matplotlib.pyplot as plt
from MeshUtils import loadVTK

# TODO
# verify derivates numerically
# spot check matrix elements (and rhs)
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
c = 3e8
w = 3e9
k = w / c

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
	R = np.linalg.norm( x - y )
	return np.exp( 1j * k * R ) / ( 4 * np.pi * R + 1e-15 )

def gradxG( x,y ):
	dx = x-y
	R = np.linalg.norm( dx )
	return dx * ( 1j * k * R - 1 ) * np.exp( 1j * k * R ) / ( 4 * np.pi * R**3 + 1e-15 )

def gradygradxG( x,y ):
	dx = x-y
	R = np.linalg.norm( dx )
	return ( np.outer( dx,dx ) * (R**2*k**2 + 3j*R*k - 3) + (R**2 - 1j*R**3*k)*np.eye(3) ) * np.exp( 1j * k * R ) / ( 4 * np.pi * R**5 + 1e-15 )

def DoubleIntegrate( x,wx, y,wy, f ): # TODO may have to manually evaluate at all point combinations
	return wx @ f( x,y ) @ wy


A = np.zeros((2*nFacets,2*nFacets),dtype=np.complex128)
b = np.zeros((2*nFacets,len(scenarios)),dtype=np.complex128)
for i in range(2*nFacets):
	for j in range(2*nFacets):
		facet_i = i % nFacets
		facet_j = j % nFacets
		# integral of  jwu G(x_i,x_j) vi^Tv_j + vi^T dxi dxj G(x_i,x_j) v_j / jwe
		I1 = 1j * w * mu * (bases[i] @ bases[j]) * DoubleIntegrate( pts[facet_i],wts[facet_i], pts[facet_j],wts[facet_j], G )
		if facet_i == facet_j : # Singularity!
			# integral of -vi^T gradG vj^T nHat
			# TODO
			I2 = 
		else:
			I2 = DoubleIntegrate( pts[facet_i],wts[facet_i], pts[facet_j],wts[facet_j], lambda x,y: bases[i] @ gradygradxG(x,y) @ bases[j] )
		I2 /= 1j * w * ep
		A[i,j] = I1 + I2
		
	for j in range(len(scenarios)):
		# integral of vi^T scenarios[j].excitation( xi )
		b[i,j] = wts[facet_i].dot( scenarios[j].excitation(pts[facet_i]).dot(bases[facet_i]) )

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

