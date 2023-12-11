
import numpy as np
from re import search

meshFileName = "SphereMesh.vtk"

f=open(meshFileName, "r")

# Get nPts
m=None
while m is None:
	m = search( "POINTS ([0-9]+)", f.readline() )
nPts = int(m[1])

# Get Pts
vertices = np.zeros((nPts,3))
for i in range(nPts):
	vertices[i,:] = [ float(x) for x in f.readline().split() ]

# Get nCells
m=None
while m is None:
	m = search( "CELLS ([0-9]+)", f.readline() )
nCells = int(m[1])

cells=[]
for i in range(nCells):
	tmp = [int(x) for x in f.readline().split()]
	if tmp[0] == 3:
		cells.append(tmp[1:])

nCells=len(cells)
cells = np.row_stack(cells)

f.close()
