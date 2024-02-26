
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
		return self.polV * np.exp( 1j * k * xyz @ self.prop )

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


dx = (pts.reshape((-1,3,1))-pts.reshape((-1,3,1)).T).reshape(pts.shape+(3,)+pts.shape[::-1]).transpose([0,1,4,3,2])
R = np.linalg.norm( dx, axis=-1 )
I_G = np.einsum('ijkl,ik,jl->ij', np.exp( 1j * k * R ) / ( 4 * np.pi * R + 1e-15 ), wts,wts )

I_gG = ... # TODO

# TODO zero out diagonal, fill in functions of R
I_ggG = np.einsum('ijklm,ijkln,ijkl,ik,jl->ijmn', dx,dx, f(...), wts,wts )
I_ggG += np.einsum('ijklm,ijkln,ijkl,mn,ik,jl->ijmn', dx,dx, g(...), np.eye(3), wts,wts )
I_ggG /= 1j * w * epsilon_0

A = 1j * w * mu_0 * bases @ bases.T
b = np.zeros((2*nFacets,len(scenarios)),dtype=np.complex128)
for i in range(2):
	for j in range(2):
		A[i*nFacets:(i+1)*nFacets,j*nFacets:(j+1)*nFacets] *= I_G
		A[i*nFacets:(i+1)*nFacets,j*nFacets:(j+1)*nFacets] += np.einsum('ijkl,ik,jl->ij' I_ggG, bases[i*nFacets:(i+1)*nFacets], bases[j*nFacets:(j+1)*nFacets])

	for j in range(len(scenarios)):
		b[i*nFacets:(i+1)*nFacets,j] = np.einsum( 'ijk,ij,ik->i', scenarios[j].excitation(pts), wts, bases[i*nFacets:(i+1)*nFacets] )

# Factor and solve
sols = np.linalg.solve( A,b )


# post process
mags = [[np.linalg.norm( np.exp( 1j * k * srcPts.dot( obs ) ) @ sols[:,3*j:3*(j+1)] ) for obs in scene.observations] for j,scene in enumerate(scenarios)]

# Plot results
#ax = plt.figure().add_subplot(projection='3d')
#ax.quiver( testPts[:,0],testPts[:,1],testPts[:,2], np.real(b[:,0]),np.real(b[:,1]),np.real(b[:,2]), length=0.1, normalize=True)
#plt.show()

plt.plot(np.linspace(0,2*np.pi,40),mags[0])
plt.show()

