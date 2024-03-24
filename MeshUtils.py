
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
	
	return vertices, facets, areas, normals, v1, v2
