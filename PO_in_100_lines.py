
import numpy as np
import matplotlib.pyplot as plt
from MeshUtils import loadVTK

# TODO
# Plane, Wedge, No collisions
# 	power?
# currents
# initial rays

# Load mesh
meshFileName = "SphereMesh.vtk"
vertices, facets = loadVTK( meshFileName )

# Compute Normals
facetNormals = np.cross( vertices[facets[:,1],:] - vertices[facets[:,0],:], vertices[facets[:,2],:] - vertices[facets[:,0],:] )
facetAreas = np.linalg.norm(facetNormals,axis=1)
for vec,norm in zip(facetNormals,facetAreas):
	vec /= norm


# Set up Scenario(s)
c = 3e8
w = 3e9
k = w / c

class PlaneWaveScenario:
	def __init__(self, propDir, polVec, obs=None):
		self.prop = propDir.copy()
		self.polV = prolVec.copy()
		
		if obs is None:
			self.observations = [-self.prop]
		else:
			self.observations = obs.copy()
	
	def excitation(self, xyz ):
		return self.polV * np.exp( 1j * k * self.prop.dot( xyz ) )

scenarios = [PlaneWaveScenario( np.array([1,0,0]), np.array([0,0,1]), [ np.array([np.cos(th),np.sin(th),0]) for th in np.linspace(0,np.pi,40) ] )]
scenarios.extend([ PlaneWaveScenario( np.array([np.cos(th),np.sin(th),0]), np.array([0,0,1]) ) for th in np.linspace(0,np.pi,5) ])


# TODO get initial rays
totRays = np.zeros( len(scenarios), dtype=np.int64 )
rayScenarioIDs = ...
rayParentIDs = ...
rayIDs = ...
rayOrigins = ...
rayDirs = ...

while nRays > 0 :
	# trace rays in queue to get collisions
	rayEndFaces = mesh.ray.intersects_first( ray_origins=rayOrigins, ray_directions=rayDirs )
	
	# categorize collisions based on (face, edge, wrong side, material, etc)
	collisionHandlerInds = getCollisionHandlers( rayEndFaces )
	newOrder = np.argsort( collisionHandlerInds )
	for arr in [rayScenarioIDs,rayParentIDs,rayIDs,rayOrigins,rayDirs,rayEndFaces]:
		arr[:] = arr[newOrder]
	# get ranges for each handler
	handlerRanges = np.zeros( nHandlers+1, dtype=np.int32 )
	handlerRanges[1:] = np.cumsum( np.bincount(collisionHandlerInds) )
	
	# Handle collisions
	for ind in range(nHandlers):
		newRays, currents = collisionHandler[ind]( handlerRanges[ind],handlerRanges[ind+1], ... )
		# store currents
		# add new rays to ray queue

	# save ray history
	# refresh ray queue (including set rayIDs)

# for each current
	# for each observation in scenarios[ current.scenario ]
		# observation += contribution from current

