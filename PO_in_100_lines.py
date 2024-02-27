
import numpy as np
import matplotlib.pyplot as plt
from MeshUtils import loadVTK

# Load mesh
meshFileName = "SphereMesh.vtk"

vertices, facets = loadVTK( meshFileName )
nVertices = len(vertices)
nFacets = len(facets)

# Compute Normals
v1 = vertices[facets[:,1],:] - vertices[facets[:,0],:]
v2 = vertices[facets[:,2],:] - vertices[facets[:,0],:]
facetNormals = np.cross(v1, v2)
facetAreas = np.linalg.norm(facetNormals,axis=1)
for faceInd in range(nFacets):
	facetNormals[faceInd,:] /= facetAreas[faceInd]


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


# get initial rays
totRays = np.zeros( len(scenarios), dtype=np.int64 )
rayScenarioIDs = ...
rayParentIDs = ...
rayIDs = ...
rayOrigins = ...
rayDirs = ...
rayEndFaces = None

# while ray queue is nonempty

	# trace rays in queue to get collisions
	rayEndFaces = mesh.ray.intersects_first( ray_origins=rayOrigins, ray_directions=rayDirs )
	
	# save ray history
	
	raysThatHit = np.where( rayEndFaces >= 0 )[0]
	
	# categorize collisions (face, edge, wrong side, material, etc)
	# clear ray queue
	
	# Handle collision
		# compute current, and store with correct scenario
		# produce new rays and add to ray queue

# for each scenario
	# for each observation
		# for each current associated with this scenario
			# observation += contribution from current
# ALT
# for each current
	# for each observation in scenarios[ current.scenario ]
		# observation += contribution from current

