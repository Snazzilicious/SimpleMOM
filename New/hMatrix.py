
# TODO
# Constructing hMatrices
#    Error checking
#    reblockify
#    deep copy
#    support size 0
# GEMM
#    check for same basis
#    revise order of low rank products
# use correct type of exceptions
# Error checking
# Testing
	# TRSM Uppers are not as accurate

import numpy as np


'''Block boundary routines
'''
def num_blocks( min_block_size, length ):
	return max( 1, length // min_block_size )

def block_size( min_block_size, length, block_index ):
	if length < min_block_size :
		return length
	remainder = length % min_block_size
	return min_block_size + ( block_index < remainder )

def block_begin( min_block_size, length, block_index ):
	remainder = length % min_block_size
	return min_block_size * block_index + min( remainder, block_index )

def block_index( min_block_size, length, position ):
	remainder = length % min_block_size
	position_part1 = min( position, remainder * (min_block_size+1) )
	position_part2 = position - position_part1
	return ( position_part1 // (min_block_size+1) ) + ( position_part2 // min_block_size )

def on_block_boundary( min_block_size, length, position ):
	if position == length :
		return True
	index = block_index( min_block_size, length, position )
	return position == block_begin( min_block_size, length, index )


def get_limits( s, n ):
	
	if not isinstance( s, slice ):
		raise TypeError("hMatrices can only be indexed with slices (':')")
	
	if s.start is not None and not -int(n) <= s.start < n :
		raise ValueError(f"Start index {s.start} is out of bounds for axis with size {n}")
	if s.stop is not None and not -int(n) < s.stop <= n :
		raise ValueError(f"End index {s.stop} is out of bounds for axis with size {n}")
	
	begin, end, stride = s.indices(n)
	
	if stride != 1 :
		raise ValueError("hMatrices do not support strided slices")
	if begin >= end :
		raise ValueError(f"Slice limits out of order {begin}:{end}")
	
	return begin,end
	
def parse_slice_key( key, nrows, ncols ):
	if len(key) != 2:
		raise ValueError(f"Incorrect number of indices for matrix: {len(key)} were indexed.")
	
	row_begin, row_end = get_limits( key[0], nrows )
	col_begin, col_end = get_limits( key[1], ncols )
	
	return row_begin, row_end, col_begin, col_end



class TreeIterator :
	def __init__( self, parent_matrix, node_index ):
		self._parent_matrix = parent_matrix
		self._node_index = node_index
	
	@property
	def shape( self ):
		return self._parent_matrix.shape( self )
	
	@property
	def dtype( self ):
		return self._parent_matrix.dtype
	
	def min_block_size( self ):
		return self._parent_matrix.min_block_size()
	
	def __getitem__( self, key ):
		(nr,nc) = self.shape
		row_begin, row_end, col_begin, col_end = parse_slice_key( key, nr, nc )
		return self.slice( row_begin, row_end, col_begin, col_end )

	# pass self to corresponding method of _parent_matrix
	def __getattr__( self, attr ):
		return lambda *args, **kwargs : getattr( self._parent_matrix, attr )( self, *args, **kwargs )
		
		

class Slice :
	def __init__( self, parent_matrix, row_begin, row_end, col_begin, col_end ):
		self._parent_node = parent_matrix
		self._row_begin = row_begin
		self._row_end = row_end
		self._col_begin = col_begin
		self._col_end = col_end
	
	def get_parent( self ):
		return self._parent_node
	
	@property
	def shape( self ):
		return ( self._row_end-self._row_begin, self._col_end-self._col_begin )
	
	@property
	def dtype( self ):
		return self._parent_node.dtype
	
	def rank( self ):
		return self._parent_node.rank()
	
	def min_block_size( self ):
		return self._parent_node.min_block_size()
	
	def block_type( self ):
		return self._parent_node.block_type()
	
	def get_row_bounds( self ):
		parent_bounds = self._parent_node.get_row_bounds()
		bounds = np.array( [self._row_begin] + [i for i in parent_bounds if self._row_begin < i < self._row_end] + [self._row_end], dtype=np.uint64 )
		return bounds - self._row_begin
	
	def get_col_bounds( self ):
		parent_bounds = self._parent_node.get_col_bounds()
		bounds = np.array( [self._col_begin] + [i for i in parent_bounds if self._col_begin < i < self._col_end] + [self._col_end], dtype=np.uint64 )
		return bounds - self._col_begin
	
	def slice( self, row_begin, row_end, col_begin, col_end ):
		row_begin += self._row_begin
		row_end += self._row_begin
		col_begin += self._col_begin
		col_end += self._col_begin
		return self._parent_node.slice( row_begin, row_end, col_begin, col_end )
	
	def __getitem__( self, key ):
		(nr,nc) = self.shape
		row_begin, row_end, col_begin, col_end = parse_slice_key( key, nr, nc )
		return self.slice( row_begin, row_end, col_begin, col_end )
	
	def get_dense_data( self ):
		D = self._parent_node.get_dense_data()
		return D[self._row_begin:self._row_end,self._col_begin:self._col_end]
	
	def get_lowrank_data( self ):
		L,R = self._parent_node.get_lowrank_data()
		return ( L[self._row_begin:self._row_end,:], R[:,self._col_begin:self._col_end] )




class hMatrix :
	def __init__( self, n_rows, n_cols, min_block_size=256, dtype=np.complex128 ):
		self._min_block_size = min_block_size
		
		self._connectivity = [None]
		self._types = ["Zero"]
		self._shapes = [ (n_rows,n_cols) ]
		self._matrices = []
		self._free_mat_inds = []
		self._dtype = dtype
	
	"""Instance methods
	"""
	def shape( self, node=None ):
		if node is None:
			return self.shape( self.root_node() )
		else:
			return self._shapes[ node._node_index ]

	@property
	def dtype( self ):
		return self._dtype
	
	def min_block_size( self ):
		return self._min_block_size
	
	def root_node( self ):
		return TreeIterator( self, 0 )
		
	def __getitem__( self, key ):
		(nr,nc) = self.shape()
		row_begin, row_end, col_begin, col_end = parse_slice_key( key, nr, nc )
		return self.slice( self.root_node(), row_begin, row_end, col_begin, col_end )
	
	def validate_iterator( self, node ):
		if not isinstance( node, TreeIterator ):
			raise TypeError("node reference is not a TreeIterator")
		if node._parent_matrix is not self:
			raise ValueError("node is not associated with this hMatrix")
	
	def rank( self, node ):
		btype = self.block_type( node )
		if btype == "Zero" :
			return 0
		elif btype == "LowRank" :
			L,_ = self.get_lowrank_data( node )
			return L.shape[1]
		else:
			raise ValueError(f"Cannot provide rank of a block of type {btype}")
	
	def block_type( self, node ):
		return self._types[ node._node_index ]
	
	"""Tree construction routines
	"""
	def partition( self, node, row_parts, col_parts ):
		self.validate_iterator( node )
		
		nrows,ncols = self.shape( node )
		block_size = self.min_block_size()
		
		if self.block_type( node ) != "Zero" :
			raise ValueError("Cannot partition non-Zero node")
		# check sorted
		if not all([ v < row_parts[i+1] for i,v in enumerate(row_parts[:-1]) ]) :
			raise ValueError("Row partitions are not strictly increasing")
		if not all([ v < col_parts[i+1] for i,v in enumerate(col_parts[:-1]) ]) :
			raise ValueError("Col partitions are not strictly increasing")
		# check in bounds
		if not all([ v >= 0 for v in row_parts ]) :
			raise ValueError("Row partition point less than 0")
		if not all([ v <= nrows for v in row_parts ]) :
			raise ValueError("Row partition point out of bounds")
		if not all([ v >= 0 for v in col_parts ]) :
			raise ValueError("Col partition point less than 0")
		if not all([ v <= ncols for v in col_parts ]) :
			raise ValueError("Col partition point out of bounds")
		# check on block boundaries
		if not all([ on_block_boundary( block_size, nrows, v ) for v in row_parts ]) :
			raise ValueError("Row partition point not on boundary")
		if not all([ on_block_boundary( block_size, ncols, v ) for v in col_parts ]) :
			raise ValueError("Col partition point not on boundary")
		
		row_parts = [0]+[i for i in row_parts if 0 < i < nrows]+[nrows]
		col_parts = [0]+[i for i in col_parts if 0 < i < ncols]+[ncols]
	
		self._connectivity[ node._node_index ] = np.zeros( ( len(row_parts)-1, len(col_parts)-1 ), dtype=np.uint64 )
		self._types[ node._node_index ] = "H"
		
		for i,(row_begin,row_end) in enumerate(zip( row_parts[:-1], row_parts[1:] )):
			for j,(col_begin,col_end) in enumerate(zip( col_parts[:-1], col_parts[1:] )):
				nrows = row_end-row_begin
				ncols = col_end-col_begin
			
				# create new zero block
				self._connectivity.append(None)
				self._types.append( "Zero" )
				self._shapes.append( (nrows,ncols) )
			
				# add edge to new child block
				self._connectivity[ node._node_index ][ i, j ] = len(self._connectivity)-1
	
	
	def insert_matrix( self, M ):
		if M.dtype != self._dtype :
			raise ValueError(f"cannot insert matrix of type {M.dtype} into hMatrix of type {self._dtype}")
		
		if len(self._free_mat_inds) > 0 :
			ind = self._free_mat_inds.pop(0)
			self._matrices[ind] = M
		else:
			ind = len(self._matrices)
			self._matrices.append(M)
		return ind
	
	def remove_matrix( self, ind ):
		self._matrices[ind] = None
		self._free_mat_inds.append(ind)
	
	def insert_dense( self, node, D=None ):
		self.validate_iterator( node )
		
		btype = self.block_type( node )
		
		if btype == "LowRank": # for converting low rank to dense
			self.remove_matrix( self._connectivity[ node._node_index ][0,0] )
			self.remove_matrix( self._connectivity[ node._node_index ][0,1] )
			self._connectivity[ node._node_index ] = None
			self._types[ node._node_index ] = "Zero"
		
		elif btype != "Zero" and btype != "Dense" :
			raise ValueError(f"Cannot insert Dense into {btype}")
		
		if (D is not None) and (D.shape != self.shape( node )) :
			raise ValueError(f"Matrix D's shape {D.shape} does not match node shape {self.shape(node)}")
		
		if D is not None :
			if self._connectivity[ node._node_index ] is None :
				mat_ind = self.insert_matrix( D )
				self._connectivity[ node._node_index ] = np.zeros( (1,1), dtype=np.uint64 )
				self._connectivity[ node._node_index ][0,0] = mat_ind
			else:
				mat_ind = self._connectivity[node._node_index][0,0]
				self._matrices[mat_ind] = D
		self._types[ node._node_index ] = "Dense"
	
	def insert_lowrank( self, node, L=None, R=None ):
		self.validate_iterator( node )
		
		btype = self.block_type( node )
		
		if btype != "Zero" and btype != "LowRank" :
			raise ValueError(f"Cannot insert LowRank into {btype}")
		
		if (L is None) ^ (R is None) :
			raise ValueError("Must insert both low rank factors")
		
		if L is not None :
			L_shape = L.shape
			R_shape = R.shape
			node_shape = self.shape( node )
			if L_shape[0] != node_shape[0] or R_shape[1] != node_shape[1] or L_shape[1] != R_shape[0]:
				raise ValueError(f"Matrix product L{L_shape} R{R_shape}, is invalid for node with shape {node_shape}")
		
		if L is not None :
			if self._connectivity[ node._node_index ] is None :
				L_ind = self.insert_matrix( L )
				R_ind = self.insert_matrix( R )
				self._connectivity[ node._node_index ] = np.zeros( (1,2), dtype=np.uint64 )
				self._connectivity[ node._node_index ][0,:] = [L_ind,R_ind]
			else:
				L_ind = self._connectivity[node._node_index][0,0]
				R_ind = self._connectivity[node._node_index][0,1]
				self._matrices[L_ind] = L
				self._matrices[R_ind] = R
		self._types[ node._node_index ] = "LowRank"
			
	
	"""Tree traversal routines
	"""
	def children_shape( self, node ):
		return self._connectivity[ node._node_index ].shape
	
	def get_child( self, node, i, j ):
		return TreeIterator( self, self._connectivity[ node._node_index ][i,j] )
	
	def get_row_bounds( self, node ):
		self.validate_iterator( node )
		
		if self.block_type( node ) != "H" :
			n,_ = self.shape( node )
			return np.array( [ 0, n ], dtype=np.uint64 )
		
		else:
			nchild = self.children_shape( node )[0]
			bounds = np.zeros( nchild+1, dtype=np.uint64 )
			for i in range(nchild):
				bounds[i+1] = bounds[i] + self.shape( self.get_child( node, i, 0 ) )[0]
			return bounds
			
	def get_col_bounds( self, node ):
		self.validate_iterator( node )
		
		if self.block_type( node ) != "H" :
			_,n = self.shape( node )
			return np.array( [ 0, n ], dtype=np.uint64 )
		
		else:
			nchild = self.children_shape( node )[1]
			bounds = np.zeros( nchild+1,  dtype=np.uint64 )
			for i in range(nchild):
				bounds[i+1] = bounds[i] + self.shape( self.get_child( node, 0, i ) )[1]
			return bounds

	def slice( self, node, row_begin, row_end, col_begin, col_end ):
		# find lowest node in tree which contains the whole slice
		fits_in_child = True
		while fits_in_child and self.block_type( node ) == "H" :
			row_bounds = self.get_row_bounds( node )
			col_bounds = self.get_col_bounds( node )
			
			# find first child intersected by slice
			row_block = np.searchsorted( row_bounds, row_begin, side='right' )-1
			col_block = np.searchsorted( col_bounds, col_begin, side='right' )-1
		
			fits_in_child = row_end <= row_bounds[row_block+1] and col_end <= col_bounds[col_block+1]
			
			if fits_in_child :
				row_begin -= row_bounds[row_block]
				row_end -= row_bounds[row_block]
				col_begin -= col_bounds[col_block]
				col_end -= col_bounds[col_block]
				
				node = self.get_child( node, row_block, col_block )
			
		return Slice( node, row_begin, row_end, col_begin, col_end )
	
	def get_dense_data( self, node ):
		if self.block_type( node ) != "Dense" :
			raise ValueError("Requested dense data from non-Dense node.")
			
		mat_ind = self._connectivity[node._node_index][0,0]
		return self._matrices[mat_ind]
	
	def get_lowrank_data( self, node ):
		if self.block_type( node ) != "LowRank" :
			raise ValueError("Requested low rank data from non-LowRank node.")
			
		L_ind = self._connectivity[node._node_index][0,0]
		R_ind = self._connectivity[node._node_index][0,1]
		return ( self._matrices[L_ind], self._matrices[R_ind] )





