
import numpy as np
import matplotlib.pyplot as plt
from re import search

# Load mesh
meshFileName = "SphereMesh.vtk"

f=open(meshFileName, "r")

# Get nVertices
m=None
while m is None:
	m = search( "POINTS ([0-9]+)", f.readline() )
nVertices = int(m[1])

# Get Vertices
vertices = np.zeros((nVertices,3))
for i in range(nVertices):
	vertices[i,:] = [ float(x) for x in f.readline().split() ]

# Get nFacets
m=None
while m is None:
	m = search( "CELLS ([0-9]+)", f.readline() )
nFacets = int(m[1])

# Get Facets
facets=[]
for i in range(nFacets):
	tmp = [int(x) for x in f.readline().split()]
	if tmp[0] == 3:
		facets.append(tmp[1:])

nFacets=len(facets)
facets = np.row_stack(facets)

f.close()

# Get Normals
bases = np.zeros((2*nFacets,3))
bases[:nFacets] = vertices[facets[:,1],:] - vertices[facets[:,0],:]
bases[nFacets:] = vertices[facets[:,2],:] - vertices[facets[:,0],:]
#facetNormals = np.cross(v1, v2)
#facetAreas = np.linalg.norm(facetNormals,axis=1)
#for faceInd in range(nFacets):
#	facetNormals[faceInd,:] /= facetAreas[faceInd]
#	v1[faceInd,:] = np.linalg.norm(v1[faceInd])
#	v2[faceInd,:] = np.linalg.norm(v2[faceInd])

# Going to use constant basis functions within facets - v1 and v2 can be vector bases
# use divergence theorem to write in terms of perimeter integral - this should entirely cancel out
# Note might not need to precompute all of the normals and such above


# Set up Scenario(s)
c = 3e8
w = 3e9
k = w / c

class PlaneWave:
	def __init__(self, propDir, polVec, obs=None):
		self.prop = propDir.copy()
		self.polV = prolVec.copy()
		
		if obs is None:
			self.observations = [-self.prop]
		else:
			self.observations = obs.copy()
	
	def excitation(self, xyz ):
		return self.polV * np.exp( 1j * k * self.prop.dot( xyz ) )

scenarios = [PlaneWave( np.array([1,0,0]), np.array([0,0,1]), [ np.array([np.cos(th),np.sin(th),0]) for th in np.linspace(0,np.pi,40) ] )]
scenarios.extend([ PlaneWave( np.array([np.cos(th),np.sin(th),0]), np.array([0,0,1]) ) for th in np.linspace(0,np.pi,5) ])


# Fill Matrix and RHS(s)
def G( x, y ):
	R = np.linalg.norm( x - y )
	return np.exp( 1j * k * R ) / ( 4 * np.pi * R + 1e-15 )

def gradgradG( x,y ):
	return 1 # TODO

def DoubleIntegrate( x,wx, y,wy, f ):
	return wx.dot( f( x,y ).dot(wy) )


A = np.zeros((2*nFacets,2*nFacets),dtype=np.complex128)
b = np.zeros((2*nFacets,len(scenarios)),dtype=np.complex128)
for i in range(2*nFacets):
	for j in range(2*nFacets):
		facet_i = i % nFacets
		facet_j = j % nFacets
		# integral of  jwu G(x_i,x_j) vi^Tv_j + vi^T dxi dxj G(x_i,x_j) v_j / jwe
		I1 = 1j * w * mu * bases[i].dot(bases[j]) * DoubleIntegrate( pts[facet_i],wts[facet_i], pts[facet_j],wts[facet_j], G )
		if facet_i == facet_j : # Singularity!
			# integral of -vi^T gradG vj^T nHat
			# TODO
			I2 = 
		else:
			I2 = DoubleIntegrate( pts[facet_i],wts[facet_i], pts[facet_j],wts[facet_j], lambda x,y: bases[i].dot(gradgradG(x,y).dot(bases[j])) )
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

