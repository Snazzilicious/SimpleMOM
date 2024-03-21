
import numpy as np
from re import search


def loadVTK( meshFileName ):
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
	
	# Facet local vector basis
	v1 = vertices[facets[:,1],:] - vertices[facets[:,0],:]
	v2 = vertices[facets[:,2],:] - vertices[facets[:,0],:]
	normals = np.cross( v1, v2 )
	
	v1_nrms = np.linalg.norm(v1,axis=1)
	areas = np.linalg.norm(normals,axis=1)
	
	for i in range(nFacets):
		normals[i] /= areas[i]
		v1[i] /= v1_nrms[i]
		
		v2[i] -= v2[i] @ v1[i]
		v2[i] /= np.linalg.norm(v2[i])
	
	return vertices, facets, areas, normals, v1, v2
