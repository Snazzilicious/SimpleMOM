

import numpy as np

from hMatrix import hMatrix


def rSVD( L, R, max_rank=None, tol=1e-5 ):

	if max_rank is None:
		max_rank = L.shape[1]
	
	QL, RL = np.linalg.qr( L, mode="reduced" )
	QR, RR = np.linalg.qr( R.T, mode="reduced" )
	
	C = RL @ RR.T
	u,s,vh = np.linalg.svd( C )
	
	new_rank = len(s)
	new_rank -= np.searchsorted( s, tol, sorter=np.arange(len(s))[::-1], side="left" )
	new_rank = max( new_rank, 1 )
	new_rank = min( new_rank, max_rank )
	
	newL = QL @ ( u[:,:new_rank] * s[:new_rank] )
	newR = vh[:new_rank,:] @ QR.T
	
	return newL, newR

def dSVD( D, max_rank=None, tol=1e-5 ):

	if max_rank is None:
		max_rank = min(D.shape)
	
	u,s,vh = np.linalg.svd( D )
	
	new_rank = len(s)
	new_rank -= np.searchsorted( s, tol, sorter=np.arange(len(s))[::-1], side="left" )
	new_rank = max( new_rank, 1 )
	new_rank = min( new_rank, max_rank )
	
	newL = u[:,:new_rank] * s[:new_rank]
	newR = vh[:new_rank,:]
	
	return newL, newR.copy()


def llaxpy( X, Y, max_rank=None, tol=1e-5 ): # TODO check for same dataset
	# y <- x + y
	if X.block_type() == "Zero" :
		return
	
	m,n = Y.shape
	y_rank = Y.rank()
	x_rank = X.rank()
	
	newL = np.zeros(( m, y_rank+x_rank ), dtype=Y.dtype )
	newR = np.zeros(( y_rank+x_rank, n ), dtype=Y.dtype )
	
	y_node = Y.get_parent()
	if Y.block_type() == "LowRank" :
		yL, yR = y_node.get_lowrank_data()
		newL[:,:y_rank] = yL[:,:]
		newR[:y_rank,:] = yR[:,:]
	
	xL, xR = X.get_lowrank_data()
	newL[ Y._row_begin:Y._row_end, y_rank: ] = xL[:,:]
	newR[ y_rank:, Y._col_begin:Y._col_end ] = xR[:,:]
	
	newL, newR = rSVD( newL, newR, max_rank, tol )
	
	k = newL.shape[1]
	if m*n <= m*k + k*n :
		D = newL @ newR
		y_node.insert_dense( D )
	else:
		y_node.insert_lowrank( newL, newR )
	

def Leaf_GEMM( alpha, A, B, C, max_rank=None, tol=1e-5 ):

	type_a = A.block_type()
	type_b = B.block_type()
	type_c = C.block_type()
	
	if (type_c == "Zero") or (type_c == "LowRank") :
		if (type_a == "Dense") and (type_b == "Dense") :
			D_a = A.get_dense_data()
			D_b = B.get_dense_data()
			
			tmp = alpha * D_a @ D_b
			
			L,R = dSVD( tmp, max_rank, tol )
		
		elif (type_a == "LowRank") and (type_b == "LowRank") :
			L_a, R_a = A.get_lowrank_data()
			L_b, R_b = B.get_lowrank_data()
			
			tmp1 = alpha * R_a @ L_b
			
			if np.prod(L_a.shape) < np.prod(R_b.shape) : # TODO
				tmp2 = L_a @ tmp1
				L = tmp2
				R = R_b
				
			else:
				tmp2 = tmp1 @ R_b
				L = L_a
				R = tmp2
			
		else:
			raise AssertionError(f"Invalid types for A and B with LowRank C: {type_a} {type_b}")
		
		wrapper = hMatrix( *(C.shape), C.min_block_size(), C.dtype )
		wrapper.insert_lowrank( wrapper.root_node(), L, R )
		llaxpy( wrapper[:,:], C, max_rank, tol )
			
	
	elif type_c == "Dense" :
		if (type_a == "Dense") and (type_b == "Dense") :
			D_a = A.get_dense_data()
			D_b = B.get_dense_data()
			D_c = C.get_dense_data()
			
			D_c += alpha * D_a @ D_b
		
		elif (type_a == "Dense") and (type_b == "LowRank") :
			D_a = A.get_dense_data()
			L_b, R_b = B.get_lowrank_data()
			D_c = C.get_dense_data()
			
			tmp = alpha * D_a @ L_b
			D_c += tmp @ R_b
		
		elif (type_a == "LowRank") and (type_b == "Dense") :
			L_a, R_a = A.get_lowrank_data()
			D_b = B.get_dense_data()
			D_c = C.get_dense_data()
			
			tmp = alpha * R_a @ D_b
			D_c += L_b @ tmp
		
		elif (type_a == "LowRank") and (type_b == "LowRank") :
			L_a, R_a = A.get_lowrank_data()
			L_b, R_b = B.get_lowrank_data()
			D_c = C.get_dense_data()
			
			tmp1 = alpha * R_a @ L_b
			
			if np.prod(L_a.shape) < np.prod(R_b.shape) : # TODO
				tmp2 = L_a @ tmp1
				D_c += tmp2 @ R_b
			else:
				tmp2 = tmp1 @ R_b
				D_c += L_a @ tmp2
				
		else:
			raise AssertionError(f"Invalid types for A and B with Dense C: {type_a} {type_b}")
		
	else:
		raise AssertionError(f"Invalid type for C matrix {type_c}")
				
			
			


def LR_to_LR_GEMM( alpha, A, B, C, max_rank=None, tol=1e-5 ):

	min_block_size = A.min_block_size()

	if A.block_type() == "LowRank" :
		
		r = A.rank()
		k = A.shape[1]
		m,n = C.shape
		
		new_a = hMatrix( r, k, min_block_size )
		a_left, a_right = A.get_lowrank_data()
		new_a.insert_dense( new_a.root_node(), a_right )
		
		new_right = hMatrix( r, n, min_block_size )
		new_right.insert_dense( new_right.root_node(), np.zeros( (r,n), dtype=C.dtype ) )
		
		hMatrixGEMM( alpha, new_a, B, new_right, max_rank, tol )
		
		lr_wrapper = hMatrix( m, n, min_block_size )
		R = new_right.root_node().get_dense_data()
		lr_wrapper.insert_lowrank( lr_wrapper.root_node(), a_left, R )
	
	else:
		
		k = B.shape[0]
		r = B.rank()
		m,n = C.shape
	
		new_b = hMatrix( k, r, min_block_size )
		b_left, b_right = B.get_lowrank_data()
		new_b.insert_dense( new_b.root_node(), b_left )
		
		new_left = hMatrix( m, r, min_block_size )
		new_left.insert_dense( new_left.root_node(), np.zeros( (m,r), dtype=C.dtype ) )
		
		hMatrixGEMM( alpha, A, new_b, new_left, max_rank, tol )
		
		lr_wrapper = hMatrix( m, n, min_block_size )
		L = new_left.root_node().get_dense_data()
		lr_wrapper.insert_lowrank( lr_wrapper.root_node(), L, b_right )
	
	llaxpy( lr_wrapper[:,:], C, max_rank, tol )
		

def queue_H_GEMM( A, B, C ):
	row_bounds = np.union1d( A.get_row_bounds(), C.get_row_bounds() )
	col_bounds = np.union1d( B.get_col_bounds(), C.get_col_bounds() )
	inr_bounds = np.union1d( A.get_col_bounds(), B.get_row_bounds() )
	
	jobs=[]
	for row_begin,row_end in zip(row_bounds[:-1],row_bounds[1:]):
		for col_begin,col_end in zip(col_bounds[:-1],col_bounds[1:]):
			for inr_begin,inr_end in zip(inr_bounds[:-1],inr_bounds[1:]):
				
				a = A[ row_begin:row_end, inr_begin:inr_end ]
				b = B[ inr_begin:inr_end, col_begin:col_end ]
				c = C[ row_begin:row_end, col_begin:col_end ]
				
				jobs.append(( a, b, c ))
	
	return jobs


def hMatrixGEMM( alpha, A, B, C, max_rank=None, tol=1e-5 ):

	A = A[:,:]
	B = B[:,:]
	C = C[:,:]

	if not (A.min_block_size() == B.min_block_size() == C.min_block_size()) :
		raise ValueError("Incompatible min_block_sizes")

	if (A.shape[0] != C.shape[0]) or (B.shape[1] != C.shape[1]) or (A.shape[1] != B.shape[0]):
		raise ValueError(f"Incompatible matrix dimensions: {A.shape}, {B.shape}, {C.shape}")
	
	if not (A.dtype == B.dtype == C.dtype) :
		raise ValueError(f"Mismatch of dtypes {A.dtype} {B.dtype} {C.dtype}")
	
	if (C.shape[0] == 0) or (C.shape[1] == 0) :
		return

	
	job_stack = [ ( A, B, C ) ]
	while len(job_stack) > 0 :

		a,b,c = job_stack.pop(0)
		if (a.shape[0] != c.shape[0]) or (b.shape[1] != c.shape[1]) or (a.shape[1] != b.shape[0]):
			raise AssertionError(f"Incompatible matrix dimensions: {a.shape}, {b.shape}, {c.shape}")
		
		type_a = a.block_type()
		type_b = b.block_type()
		type_c = c.block_type()
		
		# Trivial case
		if (type_a == "Zero") or (type_b == "Zero") :
			continue
		
		# Special case of LowRank times non-LowRank assigned to LowRank
		elif (type_c == "Zero" or type_c == "LowRank") and ((type_a == "LowRank") ^ (type_b == "LowRank")) :
			LR_to_LR_GEMM( alpha, a, b, c )
		
		# Base case
		elif (type_a != "H") and (type_b != "H") and (type_c != "H"):
			Leaf_GEMM( alpha, a, b, c )

		# General case
		else :
			job_stack[0:0] = queue_H_GEMM( a, b, c )
		
		

