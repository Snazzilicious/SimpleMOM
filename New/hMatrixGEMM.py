
def rSVD( L, R, max_rank, tol=1e-5 ):
	
	QL, RL = np.linalg.qr( L, mode="reduced" )
	QR, RR = np.linalg.qr( R.T, mode="reduced" )
	
	C = RL @ RR.T
	u,s,vh = np.linalg.svd( C )
	
	new_rank = len(s)
	new_rank -= np.searchsorted( s, tol, sorter=np.arange(len(s))[::-1], side="left" )
	new_rank = max( new_rank, 1 )
	new_rank = min( new_rank, max_rank )
	
	newL = QL @ ( u[:,:new_rank] * s[:new_rank] )
	newR = vh[:new_rank,:] @ QR
	
	return newL, newR

def rSVD( D, max_rank, tol=1e-5 ):
	
	u,s,vh = np.linalg.svd( D )
	
	new_rank = len(s)
	new_rank -= np.searchsorted( s, tol, sorter=np.arange(len(s))[::-1], side="left" )
	new_rank = max( new_rank, 1 )
	new_rank = min( new_rank, max_rank )
	
	newL = u[:,:new_rank] * s[:new_rank]
	newR = vh[:new_rank,:]
	
	return newL, newR.copy()


def llaxpy( X, Y ): # TODO Handle C is zero, check for same parent_node
	# y <- x + y
	m,n = Y.shape
	y_node = Y.get_parent()
	yL, yR = y_node.get_lowrank_data()
	y_rank = yL.shape[1]

	xL, xR = X.get_lowrank_data()
	x_rank = xL.shape[1]
	
	newL = np.zeros(( m, y_rank+x_rank ))
	newL[:,:y_rank] = yL[:,:]
	newL[ Y.row_begin:Y.row_end, y_rank: ] = xL[:,:]
	
	newR = np.zeros(( y_rank+x_rank, n ))
	newR[:y_rank,:] = yR[:,:]
	newR[ y_rank:, Y.col_begin:Y.col_end ] = xR[:,:]
	
	y_node.insert_lowrank( newL, newR ) # TODO
	

def Leaf_GEMM( alpha, A, B, C ):
	# TODO
	type_a = a.block_type()
	type_b = b.block_type()
	type_c = c.block_type()


# TODO
def LR_to_LR_GEMM( alpha, A, B, C ):

	if A.block_type() == "LowRank" :
	
		new_a.dense = A.right
		new_right.dense = zeros
		
		hMatrixGEMM( alpha, new_a, B, new_right )
		
		lr_wrapper = ( A.left, new_right )
	
	else :
	
		new_b.dense = B.left
		new_left.dense = zeros
		
		hMatrixGEMM( alpha, A, new_b, new_left )
		
		lr_wrapper = ( new_left, B.right )
	
	C = zero_to_low_rank( C )
	llaxpy( lr_wrapper, C )
	rSVD( C )
	convert_to_dense( C.get_parent() )
		

def queue_H_GEMM( A, B, C ):
	row_bounds = np.union1d( A.row_bounds(), C.row_bounds() )
	col_bounds = np.union1d( B.col_bounds(), C.col_bounds() )
	inr_bounds = np.union1d( A.col_bounds(), B.row_bounds() )
	
	jobs=[]
	for row_begin,row_end in zip(row_bounds[:-1],row_bounds[1:]):
		for col_begin,col_end in zip(col_bounds[:-1],col_bounds[1:]):
			for inr_begin,inr_end in zip(inr_bounds[:-1],inr_bounds[1:]):
				
				a = A[ row_begin:row_end, inr_begin:inr_end ]
				b = B[ inr_begin:inr_end, col_begin:col_end ]
				c = C[ row_begin:row_end, col_begin:col_end ]
				
				jobs.append(( a, b, c ))
	
	return jobs


def hMatrixGEMM( alpha, A, B, C ):

	# TODO ensure min_block_sizes match

	job_stack = [ (A[:,:],B[:,:],C[:,:]) ]
	while len(job_stack) > 0 :

		a,b,c = job_stack.pop(0)
		if a.shape[0] != c.shape[0] or b.shape[1] != c.shape[1] or a.shape[1] != b.shape[0]:
			raise ValueError(f"Incompatible matrix dimensions: {a.shape}, {b.shape}, {c.shape}")
		
		type_a = a.block_type()
		type_b = b.block_type()
		type_c = c.block_type()
		
		# Trivial case
		if type_a == "Zero" or type_b == "Zero" :
			continue
		
		# Special case of LowRank times non-LowRank assigned to LowRank
		elif (type_c == "Zero" or type_c == "LowRank") and (type_a == "LowRank" ^ type_b == "LowRank") :
			LR_to_LR_GEMM( alpha, a, b, c )
		
		# Base case
		elif type_a != "H" and type_b != "H" and type_c != "H":
			Leaf_GEMM( alpha, a, b, c )

		# General case
		else :
			job_stack[0:0] = queue_H_GEMM( a, b, c )
		
		

