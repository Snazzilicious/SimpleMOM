
# TODO
# How to know when pointing to leaf node? extra hMatrix node? add "type" labels? - leaning towards both
# Set item
#    need to identify and replace leaf nodes
#    prevent multiple processes from modifying the same root
#        most important in subroutines
#        can maybe insert zeros to separate workloads
# Add and subtract operators
#    all combinations
# Matrix multiply
#    All combinations
# Tri Solves


def parse_slice_key( key, nrows, ncols ):
	# TODO validate limits, handle 'None', negative values
	row_begin = key[0][0]
	row_end = key[0][1]
	col_begin = key[1][0]
	col_end = key[1][1]
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
		while fits_in_child:
			row_begins = node.row_begins()
			col_begins = node.col_begins()
			
			row_block = np.searchsorted( row_begins, row_begin, side='right' )-1
			col_block = np.searchsorted( col_begins, col_begin, side='right' )-1
			
			fits_in_child = row_end <= row_begins[row_block+1] and col_end <= col_begins[col_block+1]
			
			if fits_in_child and type(node.blocks[row_block][col_block]) == hMatrix:
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




class ZeroMatrix:
	def __init__( self, nrows, ncols ):
		self.nrows = nrows
		self.ncols = ncols
	
	def row_begins( self ):
		return np.array([ 0, self.nrows ])
	
	def col_begins( self ):
		return np.array([ 0, self.ncols ])
	
	@property
	def shape( self ):
		return ( self.nrows, self.ncols )


class DenseMatrix:
	def __init__( self, nrows, ncols ):
		self.dat = None
		self.nrows = nrows
		self.ncols = ncols
	
	def row_begins( self ):
		return np.array([ 0, self.nrows ])
	
	def col_begins( self ):
		return np.array([ 0, self.ncols ])
	
	@property
	def shape( self ):
		return ( self.nrows, self.ncols )


class LowRankMatrix:
	def __init__( self, nrows, ncols ):
		self.left = None
		self.right = None
		self.nrows = nrows
		self.ncols = ncols
	
	def row_begins( self ):
		return np.array([ 0, self.nrows ])
	
	def col_begins( self ):
		return np.array([ 0, self.ncols ])
	
	@property
	def shape( self ):
		return ( self.nrows, self.ncols )
	
	@property
	def rank( self ):
		return left.shape[1]
	
	

'''Single Process routines
'''

def GEMM( alpha, A, B, beta, C ):
	
	job_stack = [ (A,B,C) ]
	while len(job_stack) > 0 :
	
		a,b,c = job_stack.pop(0)
		
		type_a = a.block_type()
		type_b = b.block_type()
		type_c = c.block_type()
		
		# Catch special cases
		if type_a == "ZeroMatrix" or type_b == "ZeroMatrix":
			continue
		
		elif type_c == "ZeroMatrix" or type_c == "LowRankMatrix" 
		 and type_a == "LowRankMatrix" ^ type_b == "LowRankMatrix" 
		 and type_a == "hMatrix" ^ type_b == "hMatrix":
		 
			# augment C to receive low rank product
			ab_rank = min( a.rank(), b.rank() ) # XXX don't call rank on hmatrix
			rank_c = c.rank()
			result_rank = ab_rank + c_rank
			new_left = np.zeros((c.nrows,result_rank))
			new_right = np.zeros((result_rank,c.ncols))
			if c_rank > 0 :
				# TODO Convert c to LowRankMatrix
				new_left[:,:c.rank] = c.left[:,:]
				new_right[:c.rank,:] = c.right[:,:]
			c.left = new_left
			c.right = new_right
			
			# replace operands, with just the bases involved
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
		all_have_data = a.is_leaf() and b.is_leaf() and c.is_leaf()
		if all_have_data :
			c = alpha * a @ b + beta * c
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



def L_Solve( L, B ):
	
	job_stack = [ ( "L_Solve", L, B ) ]
	while len(job_stack) > 0 :
		
		job = job_stack.pop(0)
		
		if job[0] == "L_Solve" :
			
			l = job[1]
			b = job[2]
			
			if all_have_data :
				b = l \ b
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
			GEMM( -1.0, a, b, 1.0, c )

def U_Solve( U, B ):	

def LU_Factor( A ):
	
	

def ACA_Fill( A, fill_func ):


'''Parallelified OOC Routines
'''
def OOC_GEMM
	# Create C in memory
	
	GEMM( ... )
	
	# Write C out
	# Erase any unused datasets
	
