
# TODO
# Constructing hMatrices
#    Assigning matrix data to Dense and Low rank blocks
# How to know when pointing to leaf node? extra hMatrix node? add "type" labels? - leaning towards both
# Set item
#    need to identify and replace leaf nodes
#    LR should check if has same dataset, maybe
# Matrix multiply
#    All combinations
#        LR+LR should still check for same dataset
#    rSVD
# Tri Solves

def index_to_limits( ind, n ):
	if isinstance( ind, slice ):
		if ind.step is not None:
			"hMatrices do not support stepped slices."
		begin = 0 if ind.start is None else ind.start
		end = n if ind.stop is None else ind.stop
	elif isinstance( ind, (int,np.integer) ):
		begin = ind
		end = ind+1
	else:
		"Only integers and slices (`:`) are valid indices."
	
	if not -n <= begin < n:
		f"Start index {begin} is out of bounds for axis with size {n}"
	if not -n < end <= n:
		f"End index {end} is out of bounds for axis with size {n}"
	
	begin = n-begin if begin < 0 else begin
	end = n-end if end < 0 else end
	
	return begin,end
	
def parse_slice_key( key, nrows, ncols ):
	if len(key) != 2:
		f"Incorrect number of indices for matrix: {len(key)} were indexed."
	
	row_begin, row_end = index_to_limits( key[0], nrows )
	col_begin, col_end = index_to_limits( key[1], ncols )
	return row_begin, row_end, col_begin, col_end


class hMatrix:
	def __init__( self, nrows, ncols ):
		self.blocks = [[ ZeroMatrix(nrows,ncols) ]]
		self.nrows = nrows
		self.ncols = ncols
	
	def row_begins( self ):
		return np.cumsum( [0] + [ i[0].shape[0] for i in self.blocks ] )
	
	def col_begins( self ):
		return np.cumsum( [0] + [ i.shape[1] for i in self.blocks[0] ] )
	
	@property
	def shape( self ):
		return ( self.nrows, self.ncols )
	
	def __getitem__( self, key ):
		
		row_begin, row_end, col_begin, col_end = parse_slice_key( key, self.nrows, self.ncols )
		
		# find lowest node in the tree which completely contains the slice
		node = self
		fits_in_child = True
		while fits_in_child and node.block_type() == "hMatrix":
			row_begins = node.row_begins()
			col_begins = node.col_begins()
			
			row_block = np.searchsorted( row_begins, row_begin, side='right' )-1
			col_block = np.searchsorted( col_begins, col_begin, side='right' )-1
			
			fits_in_child = row_end <= row_begins[row_block+1] and col_end <= col_begins[col_block+1]
			
			if fits_in_child :
				row_begin -= row_begins[row_block]
				row_end -= row_begins[row_block]
				col_begin -= col_begins[col_block]
				col_end -= col_begins[col_block]
				
				node = node.blocks[row_block][col_block]
				
		return hMatrixSlice( node, row_begin, row_end, col_begin, col_end )
			
	# TODO
	def __setitem__( self, key, val ):
	
		row_begin, row_end, col_begin, col_end = parse_slice_key( key, self.nrows, self.ncols )
		
		# save children pointers b/c cannot access them through slice's root ptr once changed
		# find intersection of partitions
		# create new nodes from those that must be sliced
		# recursively set the right values in the new nodes


class hMatrixSlice:
	def __init__( self, root_node, row_begin,row_end, col_begin,col_end ):
		self.root = root_node
		self.row_begin = row_begin
		self.row_end = row_end
		self.col_begin = col_begin
		self.col_end = col_end
	
	def row_begins( self ):
		begins = [i for i in self.root.row_begins() if self.row_begin < i < self.row_end ]
		return np.cumsum( [self.row_begin] + begins + [self.row_end] ) - self.row_begin
	
	def col_begins( self ):
		begins = [i for i in self.root.col_begins() if self.col_begin < i < self.col_end ]
		return np.cumsum( [self.col_begin] + begins + [self.col_end] ) - self.col_begin
	
	@property
	def shape( self ):
		return ( self.row_end-self.row_begin, self.col_end-self.col_begin )
		
	def __getitem__( self, key ):
		
		row_begin, row_end, col_begin, col_end = parse_slice_key( key, self.nrows, self.ncols )
		
		row_begin += self.row_begin
		row_end += self.row_begin
		col_begin += self.col_begin
		col_end += self.col_begin
		
		return root[ row_begin:row_end, col_begin:col_end ]
	
	def __setitem__( self, key, val ):
		
		row_begin, row_end, col_begin, col_end = parse_slice_key( key, self.nrows, self.ncols )
		
		row_begin += self.row_begin
		row_end += self.row_begin
		col_begin += self.col_begin
		col_end += self.col_begin
		
		root[ row_begin:row_end, col_begin:col_end ] = val



	
def create_zero_block( nrows, ncols ):
	return # TODO

def create_dense_block( nrows, ncols ):
	return # TODO

def create_low_rank_block( nrows, ncols, rank ):
	return # TODO



'''Single Process routines
'''

def augment_low_rank( A, new_rank ):
	tmp = create_low_rank_block( A.shape[0], A.shape[1], new_rank )
	if A.block_type() == "LowRankMatrix":
		A_rank = A.rank()
		tmp.blocks[0][:,:A_rank] = A.blocks[0][:,:]
		tmp.blocks[1][:A_rank,:] = A.blocks[1][:,:]
	A[:,:] = tmp[:,:]
	

def GEMM( alpha, A, B, C ):
	
	if alpha == 0:
		return
	
	job_stack = [ (A,B,C) ]
	while len(job_stack) > 0 :
	
		a,b,c = job_stack.pop(0)
		if a.shape[0] != c.shape[0] or b.shape[1] != c.shape[1] or a.shape[1] != b.shape[0]:
			f"Incompatible matrix dimensions: {a.shape}, {b.shape}, {c.shape}"
		
		type_a = a.block_type()
		type_b = b.block_type()
		type_c = c.block_type()
		
		# Catch special cases
		if type_a == "ZeroMatrix" or type_b == "ZeroMatrix":
			continue
		
		elif type_c == "ZeroMatrix" or type_c == "LowRankMatrix" 
		 and type_a == "LowRankMatrix" ^ type_b == "LowRankMatrix" 
		 and type_a == "hMatrix" ^ type_b == "hMatrix":
		 	# TODO subroutinify
			# augment C to receive low rank product
			ab_rank = min( a.rank(), b.rank() ) # XXX don't call rank on hmatrix
			c_rank = c.rank()
			result_rank = ab_rank + c_rank
			augment_low_rank( c, result_rank )
			
			# replace operands, with just the bases involved TODO proper wrapping
			c_wrapper = hMatrix()
			if type_a == "LowRankMatrix":
				c.left[:,c_rank:] = a.left[:,:]
				c_wrapper.d = c.right[:,c_rank:]
				
				a_wrapper = hMatrix()
				a_wrapper.d = a.right[:,:]
				a=a_wrapper[:,:]
			else:
				c.right[c_rank:,:] = b.right[:,:]
				c_wrapper.d = c.left[c_rank:,:]
				
				b_wrapper = hMatrix()
				b_wrapper.d = b.left[:,:]
				b=b_wrapper[:,:]
			c=c_wrapper[:,:]
			
			job_stack.insert( 0, (c,"rSVD") ) # TODO handle this
		
		# Standard cases
		if type_a != "hMatrix" and type_b != "hMatrix" and type_c != "hMatrix" :
			c = alpha * a @ b + c # TODO
		else:
			row_begins = np.union1d( a.row_begins(), c.row_begins() )
			col_begins = np.union1d( b.col_begins(), c.col_begins() )
			inner_begins = np.union1d( a.col_begins(), b.row_begins() )
			
			new_jobs=[]
			for row_begin,row_end in zip(row_begins[:-1],row_begins[1:]):
				for col_begin,col_end in zip(col_begins[:-1],col_begins[1:]):
					for inner_begin,inner_end in zip(inner_begins[:-1],inner_begins[1:]):
						
						new_a = a[row_begin:row_end,inner_begin:inner_end]
						new_b = b[inner_begin:inner_end,col_begin:col_end]
						new_c = c[row_begin:row_end,col_begin:col_end]
						
						# Partition c if needed
						if new_c.node == c.node :
							c[row_begin:row_end,col_begin:col_end] = new_c
						new_c = c[row_begin:row_end,col_begin:col_end]

						new_jobs.append(( new_a, new_b, new_c ))
			
			job_stack[0:0] = new_jobs


def Leaf_L_Solve( L, B, diag ):
	if diag != "N":
		return
	if B.block_type() == "ZeroMatrix":
		return
	if L.block_type() != "DenseMatrix":
		f"Cannot solve with nondense block type"
	
	if B.block_type() == "DenseMatrix"
		rhs = B.blocks
	elif B.block_type() == "LowRankMatrix":
		rhs = B.blocks[0]
	
	M = L.blocks[:,:]
	rhs[:,:] = np.linalg.solve( L, rhs[:,:] )
	

def L_Solve( L, B, diag="N" ):
	
	job_stack = [ ( "L_Solve", L, B ) ]
	while len(job_stack) > 0 :
		
		job = job_stack.pop(0)
		
		if job[0] == "L_Solve" :
			
			l = job[1]
			b = job[2]
			# TODO catch the LowRank-hMatrix case ... maybe it doesn't need to be
			if b.block_type() == "ZeroMatrix":
				continue
			
			if l.block_type() != "hMatrix" and b.block_type() != "hMatrix" :
				
				b = l \ b # TODO
			else:
				diag_begins = np.union1d( l.col_begins(), b.row_begins() )
				col_begins = b.col_begins()
		
				new_jobs=[]
				for diag_begin,diag_end in zip(diag_begins[:-1],diag_begins[1:]):
					for col_begin,col_end in zip(col_begins[:-1],col_begins[1:]):
						new_l = l[diag_begin:diag_end,diag_begin:diag_end]
						new_b = b[diag_begin:diag_end,col_begin:col_end]
						
						new_jobs.append(( "L_Solve", new_l, new_b ))
					
					new_a = l[diag_end:,diag_begin:diag_end]
					new_b = b[diag_begin:diag_end,:]
					new_c = b[diag_end:,:]
					
					new_jobs.append(( "Schur", new_a, new_b, new_c ))
				
				job_stack[0:0] = new_jobs
		
		else:
			a = job[1]
			b = job[2]
			c = job[3]
			GEMM( -1.0, a, b, c )


def U_Solve( U, B ): # XXX need right and left solve

def LU_Factor( A ):
	
	

def ACA_Fill( A, fill_func ):


'''Parallelified OOC Routines
'''
def OOC_GEMM
	# Create C in memory
	
	GEMM( ... )
	
	# Write C out
	# Erase any unused datasets
	
