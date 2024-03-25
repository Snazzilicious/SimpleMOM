
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import mu_0,epsilon_0,speed_of_light
from scipy.integrate import quad_vec
from numba import njit
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
nFacets = len(facets)
# Facet local vector basis
bases = np.zeros((2*nFacets,3))
bases[:nFacets,:] = vertices[facets[:,1],:] - vertices[facets[:,0],:]
bases[nFacets:,:] = vertices[facets[:,2],:] - vertices[facets[:,0],:]


# Set up Excitations
w = 3e9
k = w / speed_of_light

class PlaneWave:
	def __init__(self, propDir, polVec):
		self.prop = propDir.copy()
		self.polV = polVec.copy()
	
	def excitation(self, xyz ):
		return self.polV * np.exp( -1j * k * self.prop.dot( xyz ) )

excitations = [ PlaneWave( np.array([np.cos(th),np.sin(th),0]), np.array([0,0,1]) ) for th in np.linspace(0,np.pi,5) ]

# Integration over triangles
# triangle defined by v1 and v2 originating at p0
def integrateOverTri( f, tri_i ):
	p0 = vertices[facets[tri_i,0]]
	v1 = bases[tri_i]
	v2 = bases[tri_i+nFacets]
	return quad_vec( lambda h1: quad_vec( lambda h2: f( p0 + h1*v1 + h2*v2 ), 0, 1-h1 )[0], 0, 1 )[0] * np.linalg.norm( np.cross(v1,v2) )

def integrateOver2Tris( f, tri_i, tri_j ):
	return integrateOverTri( lambda x: integrateOverTri( lambda y: f( x,y ), tri_j ), tri_i )

# integrands involving the Green's function
def integrand_G_far( x, rHat, basis_j ):
	return 1j * w * mu_0 * np.exp( -1j * k * x @ rHat ) * ( bases[basis_j] - (bases[basis_j] @ rHat)*rHat )

def integrand_G_near( x,y, test_i,basis_j ):
	dx = x-y
	R = np.linalg.norm( dx )
	I = (bases[test_i] @ bases[basis_j])
	return 1j * w * mu_0 * I * np.exp( -1j * k * R ) / ( 4 * np.pi * R + 1e-9 )

def integrand_dxG( x,y, normal, test_i,basis_j ):
	dx = x-y
	R = np.linalg.norm( dx )
	I = -( bases[test_i] @ dx ) * ( bases[basis_j] @ normal ) * ( 1j * k * R + 1 )
	return ( I / (1j * w * epsilon_0) ) * ( np.exp( -1j * k * R ) / ( 4 * np.pi * R**3 + 1e-9 ) )

def integrand_dydxG( x,y, test_i,basis_j ):
	dx = x-y
	R = np.linalg.norm( dx )
	I = bases[test_i] @ ( np.outer(dx,dx) * (R**2*k**2 - 3j*R*k - 3) + (R**2 + 1j*R**3*k)*np.eye(3) ) @ bases[basis_j]
	return ( I / (1j * w * epsilon_0) ) * ( np.exp( -1j * k * R ) / ( 4 * np.pi * R**5 + 1e-9 ) )

# Fill Matrix and RHSs
A = np.zeros( (2*nFacets,2*nFacets), dtype=np.complex128 )
b = np.zeros( 2*nFacets, dtype=np.complex128 )
for i in range(2*nFacets):
	for j in range(2*nFacets):		
		face_i = i % nFacets
		face_j = j % nFacets
		
		A[i,j] = integrateOver2Tris( lambda x,y: integrand_G_near( x,y, i,j ), face_i, face_j )
		
		if i == j :
			A[i,j] += 0 # integrate1D( integrand_dxG(x,y, i,j), pts[face_i], bnd_pts[face_j] ) # TODO
		else:
			A[i,j] += integrateOver2Tris( lambda x,y: integrand_dydxG( x,y, i,j ), face_i, face_j )
	
	for j,ex in enumerate(excitations):
		b[i,j] = -integrateOverTri( lambda x: ex.excitation( x ) @ bases[i], face_i )

# Factor and solve
sols = np.linalg.solve( A,b )


# post process
for basis_j in range(2*nFacets):
	face_j = basis_j % nFacets
	farfield[obs] += integrateOverTri( lambda x: integrand_G_far(x,obs_dir,basis_j), face_j ) * sol[basis_j,exc[obs]]


# Plot results
#ax = plt.figure().add_subplot(projection='3d')
#ax.quiver( testPts[:,0],testPts[:,1],testPts[:,2], np.real(b[:,0]),np.real(b[:,1]),np.real(b[:,2]), length=0.1, normalize=True)
#plt.show()

plt.plot(np.linspace(0,2*np.pi,40),mags[0])
plt.show()

