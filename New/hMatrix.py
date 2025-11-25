
# TODO
# Constructing hMatrices
#    todense()
#    Test construction
#    Error checking
#    Vectorize block boundaries
#    reblockify
# Support size zero slices
# GEMM
#    all leaf combinations
#        DDD, DDL, DLD, DLL*, LDD, LDL*, LLD, LLL
#        check for same basis
# semi-colons
# use correct type of exceptions



'''Block boundary routines
'''
def num_blocks( min_block_size, length ):
	return max( 1, length // min_block_size )

def block_size( min_block_size, length, block_index ):
	if min_block_size < length :
		return length
	remainder = length % min_block_size
	return min_block_size + ( block_index < remainder )

def block_index( min_block_size, length, position ):
	remainder = length % min_block_size
	position_part1 = min( position, remainder * (min_block_size+1) )
	position_part2 = position - position_part1
	return ( position_part1 // (min_block_size+1) ) + ( position_part2 // min_block_size )

def block_begin( min_block_size, length, block_index ):
	remainder = length % min_block_size
	return min_block_size * block_index + min( remainder, block_index ).

def on_block_boundary( min_block_size, length, position ):
	if position == length :
		return True
	index = block_index( min_block_size, length, position )
	return position == block_begin( min_block_size, length, index )


def get_limits( s, n ):
	
	if not isinstance( s, slice ):
		raise TypeError("hMatrices can only be indexed with slices (':')")
	
	if s.start is not None and not -n <= s.start < n :
		raise ValueError(f"Start index {s.start} is out of bounds for axis with size {n}")
	if s.end is not None and not -n < s.end <= n :
		raise ValueError(f"End index {s.end} is out of bounds for axis with size {n}")
	
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
	
	def min_block_size( self ):
		return self._parent_matrix.min_block_size()
	
	def __getitem__( self, key ):
		(nr,nc) = self.shape
		row_begin, row_end, col_begin, col_end = parse_slice_key( key, nr, nc )
		return self.slice( row_begin, row_end, col_begin, col_end )

	# pass self to corresponding method of _parent_matrix
	def __getattr__( self, attr ):
		return lambda *args, **kwargs : getattr( self.s, attr )( self, *args, **kwargs )
		
		

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
	
	def rank( self ):
		return self._parent_node.rank()
	
	def min_block_size( self ):
		return self._parent_node.min_block_size()
	
	def block_type( self ):
		return self._parent_node.block_type()
	
	def get_row_bounds( self ):
		parent_bounds = self._parent_node.get_row_bounds()
		bounds = np.array( [self._row_begin] + [i for i in parent_bounds if self._row_begin < i < self._row_end] + [self._row_end] )
		return bounds - self._row_begin
	
	def get_col_bounds( self ):
		parent_bounds = self._parent_node.get_col_bounds()
		bounds = np.array( [self._col_begin] + [i for i in parent_bounds if self._col_begin < i < self._col_end] + [self._col_end] )
		return bounds - self._col_begin
	
	def slice( self, row_begin, row_end, col_begin, col_end ):
		row_begin += self._row_begin
		row_end += self._row_begin
		col_begin += self._col_begin
		col_end += self._col_begin
		return self._parent_node.slice( self, row_begin, row_end, col_begin, col_end )
	
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
	def __init__( self, n_rows, n_cols, min_block_size=256 ):
		self._min_block_size = min_block_size
		
		self._connectivity = [None]
		self._types = ["Zero"]
		self._shapes = [ (n_rows,n_cols) ]
		self._matrices = []
		self._free_mat_inds = []
	
	"""Instance methods
	"""
	@property
	def shape( self ):
		return self._shapes[0]
	
	def min_block_size( self ):
		return self._min_block_size
	
	def root_node( self ):
		return TreeIterator( self, 0 )
		
	def __getitem__( self, key ):
		(nr,nc) = self.shape
		row_begin, row_end, col_begin, col_end = parse_slice_key( key, nr, nc )
		return self.slice( self.root_node(), row_begin, row_end, col_begin, col_end )
	
	def shape( self, node ):
		return self._shapes[ node._node_index ]
	
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
		# TODO validate_iterator
		# TODO ensure is Zero type node
		# TODO ensure parts are sorted, in bounds, and on boundaries
		nrows,ncols = self.shape( node )
		block_size = self.min_block_size()
		
		if 0 not in row_parts :
			row_parts.insert( 0, 0 )
		if nrows not in row_parts :
			row_parts.append( nrows )
		if 0 not in col_parts :
			col_parts.insert( 0, 0 )
		if ncols not in col_parts :
			col_parts.append( ncols )
	
		connectivity[ node._node_index ] = np.array( ( len(row_parts)-1, len(col_parts)-1 ), dtype=np.uint64 )
		types[ node._node_index ] = "H"
		
		for i,(row_begin,row_end) in enumerate(zip( row_parts[:-1], row_parts[1:] )):
			for j,(col_begin,col_end) in enumerate(zip( col_parts[:-1], col_parts[1:] )):
				nrows = row_end-row_begin
				ncols = col_end-col_begin
			
				# create new zero block
				self._connectivity.append(None)
				self._types.push_back( "Zero" )
				self._shapes.append( (nrows,ncols) )
			
				# add edge to new child block
				connectivity[ node._node_index ][ i, j ] = len(connectivity)-1
	
	
	def insert_matrix( self, M ):
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
	
	def insert_dense( self, node ):
		btype = self.block_type( node )
		if btype != "Zero" and btype != "Dense" :
			raise ValueError(f"Cannot insert Dense into {btype}")
		self._types[ node._node_index ] = "Dense"
	
	def insert_lowrank( self, node ):
		btype = self.block_type( node )
		if btype != "Zero" and btype != "LowRank" :
			raise ValueError(f"Cannot insert LowRank into {btype}")
		self._types[ node._node_index ] = "LowRank"
	
	def insert_dense( self, node, D ):
		if self.block_type( node ) == "LowRank": # for converting low rank to dense
			self.remove_matrix( self._connectivity[ node._node_index ][0,0] )
			self.remove_matrix( self._connectivity[ node._node_index ][0,1] )
			self._connectivity[ node._node_index ] = None
			self._types[ node._node_index ] = "Zero"
		
		if D.shape != self.shape( node ) :
			raise ValueError(f"Matrix D's shape {D.shape} does not match node shape {self.shape(node)}")
			
		self.insert_dense( node )
		
		if self._connectivity[ node._node_index ] is None :
			self._connectivity[ node._node_index ] = np.array( (1,1), dtype=np.uint64 )
			self._connectivity[ node._node_index ][0,0] = self.insert_matrix( D )
		else:
			mat_ind = self._connectivity[node._node_index][0,0]
			self._matrices[mat_ind] = D
	
	def insert_lowrank( self, node, L, R ):
		L_shape = L.shape
		R_shape = R.shape
		node_shape = self.shape( node )
		if L_shape[0] != node_shape[0] or R_shape[1] != node_shape[1] or L_shape[1] != R_shape[0]:
			raise ValueError(f"Matrix product L{L_shape} R{R_shape}, is invalid for node with shape {node_shape}")
		
		self.insert_lowrank( node )
		
		if self._connectivity[ node._node_index ] is None :
			self._connectivity[ node._node_index ] = np.array( (1,2), dtype=np.uint64 )
			self._connectivity[ node._node_index ][0,0] = self.insert_matrix( L )
			self._connectivity[ node._node_index ][0,1] = self.insert_matrix( R )
		else:
			L_ind = self._connectivity[node._node_index][0,0]
			R_ind = self._connectivity[node._node_index][0,1]
			self._matrices[L_ind] = L
			self._matrices[R_ind] = R
			
	
	"""Tree traversal routines
	"""
	def children_shape( self, node ):
		return connectivity[ node._node_index ].shape
	
	def get_child( self, node, i, j ):
		return TreeIterator( self, self.connectivity[ node._node_index ][i,j] )
	
	def get_row_bounds( self, node ):
		nchild = self.children_shape( node )[0]
		bounds = np.zeros(nchild+1)
		for i in range(nchild):
			bounds[i+1] = bounds[i] + self.shape( self.get_child( node, i, 0 ) )[0]
		return bounds
			
	def get_col_bounds( self, node ):
		nchild = self.children_shape( node )[1]
		bounds = np.zeros(nchild+1)
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
			
		i = connectivity[node._node_index][0,0]
		return self._matrices[i]
	
	def get_lowrank_data( self, node ):
		if self.block_type( node ) != "LowRank" :
			raise ValueError("Requested low rank data from non-LowRank node.")
			
		i_l = self._connectivity[node._node_index][0,0]
		i_r = self._connectivity[node._node_index][0,1]
		return ( self._matrices[i_l], self._matrices[i_r] )
			
			
	
