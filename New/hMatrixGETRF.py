


def Leaf_GETRF( piv, A ):
	a = A.get_dense_data()
	lu, swaps = scipy.linalg.lu_factor( a, overwrite_a=True, check_finite=False )
	
	a[:,:] = lu[:,:]
	
	piv[:] = np.arange( len(swaps) )
	for i in range(len(swaps)):
		piv[i], piv[swaps[i]] = piv[swaps[i]], piv[i]


def queue_H_GETRF( piv, A ):
	diag_bounds = a.get_row_bounds()
	
	jobs=[]
	for diag_begin,diag_end in zip(diag_bounds[:-1],diag_bounds[1:]):
		
		p = piv[ diag_begin:diag_end ]
		a = A[ diag_begin:diag_end, diag_begin:diag_end ]
		
		jobs.append(( "GETRF", p, a ))
		
		
		if diag_end != diag_bounds[-1] :
			b = A[ diag_begin:diag_end, diag_end: ]
			c = A[ diag_end:, diag_begin:diag_end ]
			d = A[ diag_end:, diag_end: ]
			
			jobs.append(( "TRSM_GEMM", p, a, b, c, d ))
		
	return jobs


def hMatrixGETRF( piv, A ):

	A = A[:,:]
	
	N = A.shape[0]
	if N != A.shape[1] :
		raise ValueError(f"Cannot factor non-square hMatrix {A.shape}")
	if len(piv) < N:
		raise ValueError("piv not long enough {len(piv)}")

	min_block_size = A.min_block_size()
	
	job_stack = [ ( "GETRF", piv[:N], A ) ]
	while len(job_stack) > 0 :
		
		job = job_stack.pop(0)
		
		if job[0] == "TRSM_GEMM" :
		
			p,a,b,c,d = job[1:]
			hMatrixTRSM( "L", "L", p, a, b )
			hMatrixTRSM( "R", "U", p, a, c )
			hMatrixGEMM( -1.0, c, b, d )
		
		else:
			p,a = job[1:]
			type_a = a.block_type()
			
			n = a.shape[0]
			if n != a.shape[1] :
				raise AssertionError("Non-square block on diagonal")
			if n != len(p) :
				raise AssertionError("pivot array length mismatch")
			
			# Base case
			if type_a == "Dense" :
				if n == min_block_size :
					Leaf_GETRF( p, a )
				
				else : # Block is bigger than min_block_size and must be subdivided
					wrapper = wrap_partition_dense( n, min_block_size, a )
					hMatrixGETRF( p, a )
						
			# Invalid case
			elif type_a == "Zero" or type_a == "LowRank" :
				raise ValueError("Rank-deficient diagonal block in GETRF")
			
			# General case
			else:
				job_stack[0:0] = queue_H_GETRF( p, a )





