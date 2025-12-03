
import numpy as np
from scipy.linalg import solve_triangular

from . import hMatrix, num_blocks, block_begin
from .hMatrixGEMM import hMatrixGEMM


def Leaf_TRSM( side, uplo, piv, A, B ):
	
	lu = A.get_dense_data()
	b_type = B.block_type()
	if b_type == "Dense" :
		b = B.get_dense_data()
	elif b_type == "LowRank" :
		if side == "L" :
			b,_ = B.get_lowrank_data()
		else:
			_,b = B.get_lowrank_data()
	elif b_type == "Zero" :
		return
	else:
		raise ValueError(f"Block type of B ({b_type}) invalid")
	
	
	if side == "L" and uplo == "L" :
		x = b[piv]
		x = solve_triangular( lu, x, trans=0, lower=True, unit_diagonal=True, overwrite_b=True, check_finite=False)
	
	elif side == "R" and uplo == "L" :
		x = b.T
		x = solve_triangular( lu, x, trans=1, lower=True, unit_diagonal=True, overwrite_b=True, check_finite=False)

		ipiv = np.argsort(piv)
		x = x[ipiv]
		x = x.T
	
	elif side == "L" and uplo == "U" :
		x = solve_triangular( lu, b, trans=0, lower=False, unit_diagonal=False, overwrite_b=True, check_finite=False)
	
	elif side == "R" and uplo == "U" :
		x = b.T
		x = solve_triangular( lu, x, trans=1, lower=False, unit_diagonal=False, overwrite_b=True, check_finite=False)
		x = x.T
	
	else:
		raise ValueError(f"Invalid side ({side}) or uplo ({uplo})")
	
	b[:] = x[:]


def wrap_partition_dense( n, min_block_size, D ):
	n_blocks = num_blocks( min_block_size, n )
	
	wrapper = hMatrix( n, n, min_block_size )
	root = wrapper.root_node()
	block_bounds = [ block_begin( min_block_size, n, i ) for i in range(n_blocks+1) ]
	root.partition( block_bounds, block_bounds )
	
	for i in range(n_blocks):
		for j in range(n_blocks):
			row_begin = block_bounds[i]
			row_end = block_bounds[i+1]
			col_begin = block_bounds[j]
			col_end = block_bounds[j+1]
			
			block = root.get_child( i,j )
			d = D[row_begin:row_end,col_begin:col_end].get_dense_data()
			block.insert_dense( d )
	
	return wrapper


def queue_H_TRSM( side, uplo, piv, A, B ):
	
	jobs=[]
	if side == "L" and uplo == "L" :
		diag_bounds = np.union1d( A.get_col_bounds(), B.get_row_bounds() )
		for diag_begin,diag_end in zip(diag_bounds[:-1],diag_bounds[1:]):
			
			lu = A[diag_begin:diag_end,diag_begin:diag_end]
			p = piv[diag_begin:diag_end]
			
			col_bounds = B.get_col_bounds()
			for col_begin,col_end in zip(col_bounds[:-1],col_bounds[1:]):
				
				b = B[diag_begin:diag_end,col_begin:col_end]
				jobs.append(( "TRSM", p, lu, b ))
			
			if diag_end != diag_bounds[-1] :
				a = A[diag_end:,diag_begin:diag_end]
				b = B[diag_begin:diag_end,:]
				c = B[diag_end:,:]
				jobs.append(( "GEMM", a, b, c ))
	
	elif side == "R" and uplo == "U" :
		diag_bounds = np.union1d( A.get_row_bounds(), B.get_col_bounds() )
		for diag_begin,diag_end in zip(diag_bounds[:-1],diag_bounds[1:]):
			
			lu = A[diag_begin:diag_end,diag_begin:diag_end]
			p = piv[diag_begin:diag_end]
		
			row_bounds = B.get_row_bounds()
			for row_begin,row_end in zip(row_bounds[:-1],row_bounds[1:]):
				
				b = B[row_begin:row_end,diag_begin:diag_end]
				jobs.append(( "TRSM", p, lu, b ))
			
			if diag_end != diag_bounds[-1] :
				a = B[:,diag_begin:diag_end]
				b = A[diag_begin:diag_end,diag_end:]
				c = B[:,diag_end:]
				jobs.append(( "GEMM", a, b, c ))
	
	elif side == "L" and uplo == "U" :
		diag_bounds = np.union1d( A.get_col_bounds(), B.get_row_bounds() )[::-1] # iterate in reverse order
		for diag_end,diag_begin in zip(diag_bounds[:-1],diag_bounds[1:]):
			
			lu = A[diag_begin:diag_end,diag_begin:diag_end]
			p = piv[diag_begin:diag_end]
			
			if diag_end != diag_bounds[0] :
				a = A[ diag_begin:diag_end, diag_end: ]
				b = B[ diag_end:, : ]
				c = B[ diag_begin:diag_end, : ]
				jobs.append(( "GEMM", a, b, c ))
			
			col_bounds = B.get_col_bounds()
			for col_begin,col_end in zip(col_bounds[:-1],col_bounds[1:]):
			
				b = B[ diag_begin:diag_end, col_begin:col_end ]
				jobs.append(( "TRSM", p, lu, b ))
	
	elif side == "R" and uplo == "L" :
		diag_bounds = np.union1d( A.get_row_bounds(), B.get_col_bounds() )[::-1] # iterate in reverse order
		for diag_end,diag_begin in zip(diag_bounds[:-1],diag_bounds[1:]):
			
			lu = A[diag_begin:diag_end,diag_begin:diag_end]
			p = piv[diag_begin:diag_end]
			
			if diag_end != diag_bounds[0] :
				a = B[ :, diag_end: ]
				b = A[ diag_end:, diag_begin:diag_end ]
				c = B[ :, diag_begin:diag_end ]
				jobs.append(( "GEMM", a, b, c ))
			
			row_bounds = B.get_row_bounds()
			for row_begin,row_end in zip(row_bounds[:-1],row_bounds[1:]):
				
				b = B[row_begin:row_end,diag_begin:diag_end]
				jobs.append(( "TRSM", p, lu, b ))
	
	else:
		raise ValueError(f"Invalid side ({side}) or uplo ({uplo})")
	
	return jobs
		


def hMatrixTRSM( side, uplo, piv, A, B ):

	A = A[:,:]
	B = B[:,:]
	
	if not side in ["L","R"] :
		raise ValueError("side must be 'L' or 'R', {side} passed in.")
	if not uplo in ["L","U"] :
		raise ValueError("side must be 'L' or 'U', {uplo} passed in.")
	
	min_block_size = A.min_block_size()
	
	if not A.min_block_size() == B.min_block_size() :
		raise ValueError("Incompatible min_block_sizes")
	
	if A.shape[0] != A.shape[1] :
		raise ValueError(f"Non-square matrix passed to TRSM {A.shape}")
	if side == "L" and A.shape[1] != B.shape[0] :
		raise ValueError(f"Incompatible matrix dimensions {A.shape} {B.shape}")
	elif side == "R" and A.shape[0] != B.shape[1] :
		raise ValueError(f"Incompatible matrix dimensions {B.shape} {A.shape}")
	
	if not A.dtype == B.dtype :
		raise ValueError(f"Mismatch of dtypes {A.dtype} {B.dtype}")
	
	if ( side == "L" and B.shape[1] == 0 ) or ( side == "R" and B.shape[0] == 0 ):
		return
	
	
	job_stack = [ ( "TRSM", piv[:], A, B ) ]
	while len(job_stack) > 0 :
		
		job = job_stack.pop(0)
		
		if job[0] == "GEMM" :
			a,b,c = job[1:]
			hMatrixGEMM( -1.0, a, b, c )
		
		else:
			p,a,b = job[1:]
			type_a = a.block_type()
			type_b = b.block_type()
			
			n = a.shape[0]
			if n != a.shape[1] :
				raise AssertionError("Non-square block on diagonal")
			if side == "L" and n != b.shape[0] :
				raise AssertionError(f"Incompatible matrix dimensions {a.shape} {b.shape}")
			if side == "R" and n != b.shape[1] :
				raise AssertionError(f"Incompatible matrix dimensions {a.shape} {b.shape}")
			
			# Trivial case
			if type_b == "Zero" :
				continue
			
			# Base case
			if type_a == "Dense" and type_b != "H":
				if num_blocks( min_block_size, n ) == 1 :
					Leaf_TRSM( side, uplo, p, a, b )
				
				else : # Block is bigger than min_block_size and must be subdivided
					wrapper = wrap_partition_dense( n, min_block_size, a )
					hMatrixTRSM( side, uplo, p, wrapper, b )
						
			# Invalid case
			elif type_a == "Zero" or type_a == "LowRank" :
				raise ValueError("Rank-deficient diagonal block in GETRF")
			
			# General case
			else:
				job_stack[0:0] = queue_H_TRSM( side, uplo, p, a, b )
















