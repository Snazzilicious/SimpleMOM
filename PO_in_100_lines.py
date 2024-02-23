
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
