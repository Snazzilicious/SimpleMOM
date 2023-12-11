
import numpy as np
import matplotlib.pyplot as plt
from re import search

# Load mesh
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

# Get Cells
cells=[]
for i in range(nCells):
	tmp = [int(x) for x in f.readline().split()]
	if tmp[0] == 3:
		cells.append(tmp[1:])

nCells=len(cells)
cells = np.row_stack(cells)

f.close()

# Get Normals
v1 = vertices[cells[:,1],:] - vertices[cells[:,0],:]
v2 = vertices[cells[:,2],:] - vertices[cells[:,0],:]
normals = np.cross(v1, v2)
for i in range(nCells):
	normals[i,:] = np.linalg.norm(normals[i])
#	v1[i,:] = np.linalg.norm(v1[i])
#	v2[i,:] = np.linalg.norm(v2[i])






testWt = np.array([1.,1.,2.])
testPts = np.row_stack([ testWt @ vertices[cells[i,:],:] for i in range(nCells) ]) / 4.0

srcWt = np.array([1.,2.,1.])
srcPts = np.row_stack([ srcWt @ vertices[cells[i,:],:] for i in range(nCells) ]) / 4.0

# Set up Excitation(s)
c = 3e8
w = 3e9
k = w / c

# ( propagation, polarization )
excitations = [ ( np.array([1,0,0]), np.array([0,0,1]) ) ]
#excitations = [ ( np.array([np.cos(th),np.sin(th),0]), np.array([0,0,1]) ) for th in np.linspace(0,2*np.pi,5) ]

# Set up observation(s)
observations = [ [ np.array([np.cos(th),np.sin(th),0]) for th in np.linspace(0,2*np.pi,40) ]  ]
#observations = [ [ -prop ] for prop,_ in excitations ]


# build Matrix and RHS(s)
A = np.zeros((nCells,nCells),dtype=np.complex128)
bs = np.zeros((nCells,3*len(excitations)),dtype=np.complex128)
for i,x in enumerate(testPts):
	for j,y in enumerate(srcPts):
		R = np.linalg.norm( x - y )
		A[i,j] = np.exp( 1j * k * R ) / ( 4 * np.pi * R )

	for j,(prop,polar) in enumerate(excitations):
		bs[i,3*j:3*(j+1)] = - polar * np.exp( 1j * k * prop.dot( x ) )
	


# Factor and solve
sols = np.linalg.solve( A,bs )

# post process
mags = [[np.linalg.norm( np.exp( 1j * k * srcPts.dot( obs ) ) @ sols[:,3*j:3*(j+1)] ) for obs in observations[j]] for j in range(len(excitations)) ]


# Plot results
#ax = plt.figure().add_subplot(projection='3d')
#ax.quiver( testPts[:,0],testPts[:,1],testPts[:,2], np.real(b[:,0]),np.real(b[:,1]),np.real(b[:,2]), length=0.1, normalize=True)
#plt.show()

plt.plot(np.linspace(0,2*np.pi,40),mags[0])
plt.show()

