"""Functions for reading, writing, and processing meshes for simulation.

Conventions:
	0-based indexing
	Arrays of vertices are ndarrays with shape ( num_vertices, ndim ), where ndim=3 in almost all circumstances
	Arrays of homogeneous elements are ndarrays with shape ( num_elements, num_vertices ), where e.g. num_vertices=3 for triangles.
"""
import numpy as np


""" VTK Files
Vertex indices are 0-indexed in the file
"""

from re import search
	
def load_legacy_vtk( meshfilename ):
	"""Loads a mesh (vertices and triangles) from a legacy VTK file
	
	Arguments
		meshfilename : string
			Name of the legacy VTK file from which to load the mesh.
	
	Returns
		vertices : ndarray of float64 with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh
		faces : ndarray of int64 with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
	"""
	f=open(meshfilename, "r")
	
	# Get n_vertices
	m=None
	while m is None:
		m = search( "POINTS ([0-9]+)", f.readline() )
	n_vertices = int(m[1])

	# Get Vertices
	vertices = np.zeros((n_vertices,3))
	for i in range(n_vertices):
		vertices[i,:] = [ float(x) for x in f.readline().split() ]

	# Get n_faces
	m=None
	while m is None:
		m = search( "CELLS ([0-9]+)", f.readline() )
	n_faces = int(m[1])

	# Get Facets
	faces=[]
	for i in range(n_faces):
		tmp = [int(x) for x in f.readline().split()]
		if tmp[0] == 3:
			faces.append(tmp[1:])

	n_faces=len(faces)
	faces = np.vstack(faces)

	f.close()
	
	return vertices, faces


def write_legacy_vtk( filename, vertices, faces, vertex_scalars={},vertex_vectors={}, face_scalars={},face_vectors={} ):
	"""Writes a mesh (vertices, triangles) and fields defined on that mesh to legacy VTK format for e.g. visualization in Paraview.
	
	Arguments
		filename : string
			Name of the VTK file to which to write the mesh. Must include the extension.
		vertices : ndarray of floats with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh
		faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
		vertex_scalars : dictionary of field name : array of field data
			Scalar fields defined at the mesh vertices. The keys of the dictionary are the names of the fields.
			The values are arrays of one float value per mesh vertex.
		vertex_vectors : dictionary of field name : ndarray (n,3) of field data
			Vector fields defined at the mesh vertices. The keys of the dictionary are the names of the fields.
			The values are ndarrays of shape ( num_vertices, 3 ) of floats defining the field.
		face_scalars : dictionary of field name : array of field data
			Scalar fields defined on the mesh faces. The keys of the dictionary are the names of the fields.
			The values are arrays of one float value per mesh vertex.
		face_vectors : dictionary of field name : ndarray (n,3) of field data
			Vector fields defined on the mesh faces. The keys of the dictionary are the names of the fields.
			The values are ndarrays of shape ( num_faces, 3 ) of floats defining the field.
	"""
	f = open(filename,'w')
	
	# header
	f.write("# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\n")
	
	# print coordinates of points
	num_verts = len(vertices)
	f.write("POINTS " + str(num_verts) + " float\n")
	for p in vertices:
		f.write( "{0} {1} {2}\n".format(p[0], p[1], p[2]) )
	
	# print coordinates of vertices
	num_faces = len(faces)
	f.write("CELLS " + str(num_faces) + " " + str(4*num_faces) + "\n")
	for c in faces:
		f.write( "3 {0} {1} {2}\n".format(c[0], c[1], c[2]) )

	# cell types
	f.write("CELL_TYPES " + str(num_faces) + "\n")
	for _ in range(num_faces):
		f.write("5\n")
	
	# actual values
	if len(vertex_scalars) > 0 or len(vertex_vectors) > 0 :
		f.write( "POINT_DATA " + str(num_verts) + "\n" )
	
		for name in vertex_scalars:
			f.write("SCALARS {0} float\nLOOKUP_TABLE default\n".format(name))
			for v in vertex_scalars[name]:
				f.write("{0}\n".format(v))
		
		for name in vertex_vectors:
			f.write("VECTORS {0} float\n".format(name))
			for v in vertex_vectors[name]:
				f.write("{0} {1} {2}\n".format(v[0],v[1],v[2]))
	
	if len(face_scalars) > 0 or len(face_vectors) > 0 :
		f.write( "CELL_DATA " + str(num_faces) + "\n" )
	
		for name in face_scalars:
			f.write("SCALARS {0} float\nLOOKUP_TABLE default\n".format(name))
			for v in face_scalars[name]:
				f.write("{0}\n".format(v))
		
		for name in face_vectors:
			f.write("VECTORS {0} float\n".format(name))
			for v in face_vectors[name]:
				f.write("{0} {1} {2}\n".format(v[0],v[1],v[2]))

	f.close()


import xml.etree.ElementTree as ET
from sys import byteorder as endianness
from base64 import b64encode

def add_data_array( parent, name, data_array, binary=True ):
	"""Helper function for write_vtk.
	"""
	darray = ET.SubElement( parent, "DataArray" )
	
	darray.set( "Name", name )
	
	if len(data_array.shape) > 1:
		darray.set( "NumberOfComponents", str(data_array.shape[-1]) )
	
	tp = str(data_array.dtype)
	if tp[0] == 'u' :
		tp = tp[:2].upper() + tp[2:]
	else:
		tp = tp[:1].upper() + tp[1:]
	darray.set( "type", tp )
	
	if binary :
		darray.set( "format", "binary" )
		num_bytes = data_array.itemsize*len(data_array.reshape(-1))
		darray.text = str(b64encode( num_bytes.to_bytes( 8, endianness ) ))[2:-1] + str(b64encode( data_array.reshape(-1).tobytes() ))[2:-1]
	else:
		darray.set( "format", "ascii" )
		darray.text = " ".join([ str(i) for i in data_array.reshape(-1) ])


def write_vtk( filename, vertices, cells, vertex_values={}, face_values={} ):
	"""Writes a mesh (vertices, triangles) and fields defined on that mesh to VTK XML format for e.g. visualization in Paraview.
	
	Arguments
		filename : string
			Name of the VTK file to which to write the mesh. Must include the extension.
		vertices : ndarray of floats with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh
		faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
		vertex_values : dictionary of field name : array of field data
			Fields defined at the mesh vertices. The keys of the dictionary are the names of the fields.
			The values are arrays field data. Scalar fields are indicated by arrays of length num_vertices of floats,
			vector fields are indicated by ndarrays of shape ( num_vertices, 3 ) of floats.
		face_values : dictionary of field name : array of field data
			Fields defined on the mesh faces. The keys of the dictionary are the names of the fields.
			The values are arrays field data. Scalar fields are indicated by arrays of length num_faces of floats,
			vector fields are indicated by ndarrays of shape ( num_faces, 3 ) of floats.
	"""
	num_verts = len(vertices)
	num_faces = len(cells)
	
	root = ET.Element( "VTKFile" )
	root.set( "type", "UnstructuredGrid" )
	root.set( "version", "1.0" )
	root.set( "byte_order", { "little" : "LittleEndian", "big" : "BigEndian" }[endianness] )
	root.set( "header_type", "UInt64" )

	ugrid = ET.SubElement( root, "UnstructuredGrid" )
	
	piece0 = ET.SubElement( ugrid, "Piece" )
	piece0.set( "NumberOfPoints", str(num_verts) )
	piece0.set( "NumberOfCells", str(num_faces) )
	
	# Mesh
	pts = ET.SubElement( piece0, "Points" )
	add_data_array( pts, "Points", vertices )
	cls = ET.SubElement( piece0, "Cells" )
	add_data_array( cls, "connectivity", cells.reshape(-1) ) # vertices indices
	add_data_array( cls, "offsets", 3*np.arange(num_faces)+3 ) # ends of cells in previous array
	add_data_array( cls, "types", 5*np.ones(num_faces, dtype=np.uint8) ) # cell types
	
	# Values
	if len(vertex_values) > 0 :
		pt_data = ET.SubElement( piece0, "PointData" )
		for name in vertex_values:
			add_data_array( pt_data, name, vertex_values[name] )
		
	if len(face_values) > 0 :
		pt_data = ET.SubElement( piece0, "CellData" )
		for name in face_values:			
			add_data_array( pt_data, name, face_values[name] )

	# Write to file
	from xml.dom.minidom import parseString
	xmlstr = parseString( ET.tostring(root, xml_declaration=True) )
	with open(filename, "w") as f:
		f.write( xmlstr.toprettyxml(indent="\t") )
	# TODO when have Python >= 3.9 do this
#	tree = ET.ElementTree( element=root )
#	ET.indent( tree, space='\t' )
#	tree.write( filename, xml_declaration=True )




""" PAT & STARS Files
Vertex indices are 1-indexed in the file
"""

def load_pat( meshfilename ):
	"""Loads a mesh (vertices, triangles, mesh groups, and group names) from a PATRAN file.
	
	Arguments
		meshfilename : string
			Name of the PATRAN file from which to load the mesh.
	
	Returns
		vertices : ndarray of float64 with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh
		faces : ndarray of int64 with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
		groups : list of list
			Members of each group of faces. groups[i] is a list of face indices belonging to group i.
		groupnames : list of string
			The names of each group.
	"""
	f=open(meshfilename,'r')
	for _ in range(4):
		f.readline()
	nextline = f.readline()

	vertices=[]
	while nextline[:2] == ' 1' :
		vertices.append([float(i) for i in f.readline().split()])
		f.readline()
		nextline = f.readline()

	vertices = np.vstack(vertices)

	# XXX These indices are 1-based XXX
	faces=[]
	while nextline[:2] == ' 2' :
		f.readline()
		faces.append([int(i) for i in f.readline().split()])
		nextline = f.readline()

	faces = np.vstack(faces)-1

	# XXX so are these XXX
	groupnames=[]
	groups=[]
	while nextline[:2] == '21' :
		n_vals = int(nextline.split()[2])
		tally=0
		groupnames.append(f.readline()[:-1])
		groups.append([])
		while tally < n_vals :
			groups[-1].extend([ int(v)-1 for i,v in enumerate(f.readline().split()) if i % 2 == 1 ])
			tally += 10
		nextline = f.readline()

	f.close()
	
	return vertices, faces, groups, groupnames


import re

def load_stars( filename ):
	"""Loads a mesh (vertices and triangles) from a .STARS file.
	
	Arguments
		meshfilename : string
			Name of the STARS file from which to load the mesh.
	
	Returns
		vertices : ndarray of float64 with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh
		faces : ndarray of int64 with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
	"""
	f=open(filename,'r')
	line = f.readline()
	
	while not "coordinates" in line:
		line = f.readline()
	n_vertices = int( line.split()[-1] )
	
	while not "triangles" in line:
		line = f.readline()
	n_faces = int( line.split()[-1] )
	
	while not "mesh IDs" in line:
		line = f.readline()
	f.readline()
	
	# get vertices
	verts_string_list = f.readline().split()
	vertices=[]
	for vert_str in verts_string_list:
		m = re.findall( r"[0\-]\.[0-9]+E[+\-][0-9][-0-9]", vert_str )
		vertices.append([ float(i) for i in m ])
	
	vertices = np.vstack(vertices)
	
	f.readline()
	
	# get faces
	faces=[]
	for _ in range(n_faces):
		faces.append([ int(i) for i in f.readline().split(',')[:3] ])
	
	faces = np.vstack(faces)-1
	
	return vertices, faces


import datetime

def write_stars( filename, vertices, faces, groups, groupnames ):
	"""Writes a mesh (vertices, triangles, and groups) to the .STARS format for simulating in Stars.
	
	Arguments
		filename : string
			Name of the file to which to write the mesh. Must include the extension.
		vertices : ndarray of floats with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh
		faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
		groups : list of list
			Members of each group of faces. groups[i] is a list of face indices belonging to group i.
		groupnames : list of string
			The names of each group. Currently unused, can be `None`.
	"""
	f=open(filename,'w')

	dt = datetime.datetime.now()
	f.write(" Title: PATRAN ACAD-12.8, {0:4d}{1:02d}{2:02d}, 152824\n".format( dt.year, dt.month, dt.day ))
	f.write(" Mesh Units:    M\n")
	f.write(" Number of coordinates:{0:15d}\n".format(len(vertices)))
	f.write(" Number of bars :{0:21d}\n".format(0))
	f.write(" Number of triangles:{0:17d}\n".format(len(faces)))
	f.write(" Number of tets:{0:22d}\n".format(0))
	f.write(" Number of mesh IDs:{0:18d}\n".format(len(groups)))
	f.write("\n")

	# printing values in STARS' ridiculous format
	vertex_str = " ".join(["{0:.6E}{1:.6E}{2:.6E}".format( vert[0],vert[1],vert[2] ) for vert in vertices])

	def swap_decimal(m):
		s=m[0]
		if float(s) > 0:
			return "0."+s[0]
		else:
			return "-."+s[1]

	vertex_str = re.sub( r"-*[1-9]\.", swap_decimal, vertex_str )

	def exp_incr(m):
		val = int(m[0][1:])+1
		if val == 0:
			return "E+00"
		else:
			return m[0][:2]+"{0:02d}".format(abs(val))

	vertex_str = re.sub( r"E[+-][0-9][0-9]", exp_incr, vertex_str )

	vertex_str = re.sub( r"0\.000000E\+01","0.0000000E+00", vertex_str )
	vertex_str = re.sub( r"-0\.","-.", vertex_str )

	f.write(vertex_str)

	f.write("\nTriangles:{0:12d}\n".format(len(faces)))

	for i,g in enumerate(groups):
		for face_ind in g:
			face = faces[face_ind]
			f.write( "{0},{1},{2},{3}\n".format( face[0]+1,face[1]+1,face[2]+1, i+1 ) )

	f.close()



"""GEMs mesh files
"""

def load_geo( filename ):
	"""Loads a mesh (vertices, triangles, and mesh groups) from a .geo file.
	
	Arguments
		meshfilename : string
			Name of the .geo file from which to load the mesh.
	
	Returns
		vertices : ndarray of float64 with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh
		faces : ndarray of int64 with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
		groups : list of list
			Members of each group of faces. groups[i] is a list of face indices belonging to group i.
	"""
	f=open(filename,"r")

	sline=f.readline().split()
	vertices=[]
	while "'Node'" in sline:
		vertices.append([float(i) for i in sline[-3:]])
		sline=f.readline().split()
	vertices=np.array(vertices)

	faces=[]
	groups=[[]]
	i=0
	currGroup=0
	while "'Triangle'" in sline:
		faces.append([int(i)-1 for i in sline[-4:-1]])
		group_i = int(sline[-1])-1
		if group_i != currGroup:
			groups.append([])
			currGroup = group_i # += 1
		groups[-1].append(i)
		i+=1
		sline=f.readline().split()
	faces=np.array(faces)

	return vertices, faces, groups


def write_geo(filename, vertices, faces, groups ):
	"""Writes a mesh (vertices, triangles, and groups) to the .geo format for simulating in GEMS.
	
	Arguments
		filename : string
			Name of the file to which to write the mesh. Must include the extension.
		vertices : ndarray of floats with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh
		faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
		groups : list of list
			Members of each group of faces. groups[i] is a list of face indices belonging to group i.
	"""
	f=open(filename,'w')

	for i,vert in enumerate(vertices):
		f.write( "'Node' {0:>21d}   {1:15f}   {2:15f}   {3:15f}\n".format( i+1, vert[0],vert[1],vert[2] ) )

	tri_i=1
	for group_i,g in enumerate(groups):
		for face_i in g:
			face = faces[face_i]
			f.write( "'Triangle' {0:>18d} 3 {1:10d} {2:10d} {3:10d} {4:10d}\n".format( tri_i, face[0]+1,face[1]+1,face[2]+1, group_i+1 ) )
			tri_i += 1
	
	f.close()



""" OBJ Files
Vertex indices are 1-indexed in the file
"""

def load_obj( meshfilename ):
	"""Loads a mesh (vertices, triangles, mesh groups, and group names) from an OBJ file.
	
	Arguments
		meshfilename : string
			Name of the OBJ file from which to load the mesh.
	
	Returns
		vertices : ndarray of float64 with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh
		faces : ndarray of int64 with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
		groups : list of list
			Members of each group of faces. groups[i] is a list of face indices belonging to group i.
		groupnames : list of string
			The names of each group.
	"""
	f = open(meshfilename,"r")
	line = f.readline()
	while line[0] != 'v' :
		line = f.readline()

	vertices=[]
	while line[:2] == 'v ' :
		vertices.append([float(i) for i in line.split()[1:]])
		line = f.readline()
	
	# XXX normals indicated by line[:2] == 'vn'

	faces=[]
	groupnames=[]
	groups=[]
	ind=0
	while len(line) > 2:
		if line[0] == 'g' :
			groupnames.append(line[2:-1])
			groups.append([])
		elif line[0] == 'f' :
			faces.append([ int(i) for i in re.findall( r"([0-9]+)(?://[0-9]+)?", line ) ])
			groups[-1].append(ind)
			ind += 1
		line = f.readline()

	f.close()

	vertices = np.vstack(vertices)
	faces = np.vstack(faces)-1
	
	return vertices, faces, groups, groupnames


def write_obj( filename, vertices, faces, groups, groupnames ):
	"""Writes a mesh (vertices, triangles, groups, and groupnames) to the OBJ format for simulating in Aurora.
	
	Arguments
		filename : string
			Name of the file to which to write the mesh. Must include the extension.
		vertices : ndarray of floats with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh
		faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
		groups : list of list
			Members of each group of faces. groups[i] is a list of face indices belonging to group i.
		groupnames : list of string
			The names of each group.
	"""
	f=open(filename,'w')

	for vert in vertices:
		f.write( "v {0} {1} {2}\n".format( vert[0],vert[1],vert[2] ) )

	for gname,g in zip(groupnames,groups):
		f.write( "g {0}\n".format(gname) )
		for face_ind in g:
			face = faces[face_ind]
			f.write( "f {0} {1} {2}\n".format( face[0]+1,face[1]+1,face[2]+1 ) )

	f.close()



""" Geometry Calculations
"""

def local_basis_areas( vertices, faces ):
	"""Computes an orthonormal vector basis for each face including the outward normal.
	Also returns the area of each face since it is a byproduct of the basis calculation.
	
	Arguments
		vertices : ndarray of floats with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh
		faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
	
	Returns
		normals : ndarray of shape ( num_faces, 3 )
			The outward, unit normal of each face in the mesh.
		v1 : ndarray of shape ( num_faces, 3 )
			First unit basis vector in the plane of the face.
			Points from the first vertex to the second.
		v2 : ndarray of shape ( num_faces, 3 )
			Second unit basis vector in the plane of the face.
			Is orthogonal to v1, but otherwise points towards the third vertex from the first.
		areas : numpy array of length num_faces of floats
			The area of each face in the mesh, in mesh units.
	"""
	# Achieves outward normal
	v1 = vertices[faces[:,1],:] - vertices[faces[:,0],:]
	v2 = vertices[faces[:,2],:] - vertices[faces[:,0],:]
	normals = np.cross( v1,v2 )
	areas = np.linalg.norm(normals,axis=1) / 2.0
	
	# normalize
	normals = ( normals.T / (2*areas) ).T
	v1 = ( v1.T / np.linalg.norm(v1,axis=1) ).T
	# and orthogonalize
	v2 = v2.T - np.sum( v1*v2, axis=1 ) * v1.T
	v2 = ( v2 / np.linalg.norm(v2,axis=0) ).T

	return normals, v1, v2, areas


def get_centroids( vertices, faces ):
	"""Computes the centroid (average of all vertices) of each face.
	
	Arguments
		vertices : ndarray of floats with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh
		faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
	
	Returns
		centroids : ndarray of shape ( num_faces, 3 )
			The centroid of each face, computed as the mean of all vertices of the face.
	"""
	centroids = vertices[faces[:,0]] + vertices[faces[:,1]] + vertices[faces[:,2]]
	return centroids / 3.0


""" Connectivity
"""
def get_nbr_map( faces ):
	"""Produces a map from each face to its neighboring faces.
	
	Arguments
		faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
	
	Returns
		nbr_map : list of lists
			nbr_map[i] is a list of indices of faces that share vertices with face i.
	"""
	# get faces attached to each vertex
	n_vertices = np.max(faces)+1
	vertex_members = [[] for _ in range(n_vertices)]
	for i,f in enumerate(faces):
		for v in f:
			vertex_members[v].append(i)
	
	# construct neighbors list
	nbr_map = [[] for _ in range(len(faces))]
	for i,f in enumerate(faces):
		for v in f:
			nbr_map[i].extend([n for n in vertex_members[v] if n != i])
		nbr_map[i] = np.unique(nbr_map[i])
	
	return nbr_map


def get_edges( faces ):
	"""Produces a map from edge to the indices of faces to which it belongs.
	
	Arguments
		faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
	
	Returns
		edges : dictionary with edges as keys and list of face indices as values.
			Maps from edges to faces including the edge. The edge is specified as a tuple of
			vertex indices defining the edge sorted least to greatest. The value is a list
			of face indices. E.G. { (min_vertex_ind,max_vertex_index) : [face1 <, face2>] }
	"""
	edges = {}
	for i,f in enumerate(faces):
		i0 = 0
		i1 = 1
		for _ in range(3):
			key = ( min( f[i0], f[i1] ), max( f[i0], f[i1] ) )
			
			edge = edges.get( key, None )
			if edge is None:
				edges[key] = [i]
			else:
				edge.append( i )
			
			i0 += 1
			i1 += 1
			i1 %= 3

	return edges


def get_edge_nbr_map( n_faces, edges ):
	"""Produces a map from each face to faces with which it shares and edge.
	
	Arguments
		n_faces : integer
			The number of faces in the mesh.
		edges : output of get_edges
			Map from edge to faces including it.
	
	Returns
		nbr_map : list of lists
			nbr_map[i] is a list of indices of faces that share edges with face i.
	"""
	nbr_map = [[] for _ in range(n_faces)]
	
	for e in edges:
		if len(edges[e]) > 1:
			i = edges[e][0]
			j = edges[e][1]
			
			nbr_map[i].append( j )
			nbr_map[j].append( i )
	
	for i in nbr_map:
		i.sort()
	
	return nbr_map





""" Finding MOM regions
"""

def get_mom_faces( nbr_map, normals, tol=0.5 ):
	"""Find faces in the mesh in regions of high curvature.
	
	Arguments
		nbr_map : list of lists
			nbr_map[i] is a list of indices of faces adjacent to face i.
		normals : ndarray of shape ( num_faces, 3 )
			The outward, unit normal of each face in the mesh.
		tol (optional) : float
			Threshold for indentifying a face as being in a region of high curvature.
			If normals[i].dot( normals[j] ) < tol for any j in nbr_map[i], then is_mom[i] is true.
	
	Returns
		is_mom : numpy array of bool
			If is_mom[i] is true, then face i is in a region of high curvature.
	"""
	n_faces = len(nbr_map)
	is_mom = np.zeros( n_faces, dtype=np.bool_ )
	
	# find where neighboring faces are not very coplanar
	for i,nbrs in enumerate(nbr_map):
		is_mom[i] = np.min(normals[nbrs,:].dot(normals[i])) < tol
	
	return is_mom


def get_mom_regions( nbr_map, normals, tol=0.5, breadth=3 ):
	"""Find connected regions in the mesh with high curvature.
	
	Arguments
		nbr_map : list of lists
			nbr_map[i] is a list of indices of faces adjacent to face i.
		normals : ndarray of shape ( num_faces, 3 )
			The outward, unit normal of each face in the mesh.
		tol (optional) : float
			Threshold for indentifying a face as being in a region of high curvature.
			If normals[i].dot( normals[j] ) < tol for any j in nbr_map[i], then curvature is considered high.
		breadth (optional) : integer
			The number of degrees of separation around a face of high curvature to include in the region.
			If breadth = 2, then the neighbors and neighbors of neighbors of each face of high curvature is
			included in the region.
	
	Returns
		mom_group_inds : numpy array of ints
			The index of the high curvature region to which each face belongs. A value of -1 indicates
			the face is not in a region of high curvature.
	"""
	is_mom = get_mom_faces( nbr_map, normals, tol )
	
	# expand mom regions
	face_stack = np.where(is_mom)[0].tolist()
	for _ in range(breadth):
		new_stack=[]
		for i in face_stack:
			for nbr in nbr_map[i]:
				if not is_mom[nbr]:
					is_mom[nbr] = True
					new_stack.append(nbr)
		face_stack = new_stack
	
	# find all connected mom regions
	mom_group_inds = -np.ones( n_faces, dtype=np.int32 )
	group_id = 0
	for i in range(n_faces):
		if is_mom[i] and mom_group_inds[i] == -1:
			nbr_stack=[i]
			while len(nbr_stack) > 0 :
				ind = nbr_stack.pop()
				if is_mom[ind] and mom_group_inds[ind] == -1:
					mom_group_inds[ind] = group_id
					nbr_stack.extend(nbr_map[ind])
			
			group_id += 1
	
	return mom_group_inds



""" Symmetry
"""

def remove_unused( vertices, faces ):
	"""Removes any unused vertices.
	
	Arguments
		vertices : ndarray of floats with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh.
		faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
	
	Returns
		new_vertices : ndarray of floats with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh,
			with unused vertices eliminated.
		new_faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
			Inidices have been adjusted to refer to the new list of vertices.
	"""
	unique_inds = np.unique(faces)
	unmap = { old_i : new_i for new_i,old_i in enumerate(unique_inds) }
	
	new_vertices = vertices[unique_inds]
	
	face_shp = faces.shape
	new_faces=faces.reshape(-1).copy()
	for i in range(len(new_faces)):
		new_faces[i] = unmap[ new_faces[i] ]
	new_faces=new_faces.reshape(face_shp)
	
	return new_vertices, new_faces
	


def mirror( plane, vertices, faces, groups, tol=1e-12 ):
	"""Mirrors a mesh about a specified plane.
	
	Arguments
		plane : "YZ", "XZ", or "XY"
			The plane about which to mirror the mesh.
		vertices : ndarray of floats with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh.
		faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
		groups : list of list
			Members of each group of faces. groups[i] is a list of face indices belonging to group i.
		tol (optional, default=1e-12) : float
			A point this close to another is considered a duplicate and is removed.
	
	Returns
		new_vertices : ndarray of floats with shape ( num_vertices, 3 )
			The xyz coordinates in space of the nodes in the mesh,
			with unused vertices eliminated.
		new_faces : ndarray of integers with shape ( num_faces, 3 )
			The indices of mesh vertices defining each triangle in the mesh.
			Inidices have been adjusted to refer to the new list of vertices.
		new_groups : list of list
			Expanded `groups` including mirrored faces.
			new_groups[i] contains each face in groups[i] and its image.
	"""
	# create mirrored points
	vertices2 = vertices.copy()
	if plane == "YZ":
		vertices2[:,0] *= -1
	elif plane == "XZ":
		vertices2[:,1] *= -1
	elif plane == "XY":
		vertices2[:,2] *= -1
	else:
		print(f"Invalid plane: {plane}")
		return None
	
	# create mirrored faces
	faces2 = faces.copy()+len(vertices)
	tmp = faces2[:,0].copy()
	faces2[:,0] = faces2[:,1]
	faces2[:,1] = tmp[:]
	
	# Expand groups to include new faces
	n_faces = len(faces)
	new_groups = [ [i for i in g] + [i+n_faces for i in g] for g in groups ]
	
	# find duplicate points
	from sklearn.neighbors import KDTree
	kdt = KDTree(vertices)
	dists,inds = kdt.query( vertices2 )
	comap = { i+len(vertices) : inds[i] for i in np.where(dists<=tol)[0] }

	# reassign where duplicate points are used
	faces2=faces2.reshape(-1)
	for i in range(len(faces2)):
		faces2[i] = comap.get( faces2[i], faces2[i] )
	faces2=faces2.reshape((-1,3))

	new_vertices = np.vstack(( vertices,vertices2 ))
	new_faces = np.vstack(( faces,faces2 ))
	new_vertices,new_faces = remove_unused( new_vertices,new_faces )
	
	return new_vertices, new_faces, new_groups

