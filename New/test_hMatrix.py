

# TODO
# TRSM
	# compare to todense tri solve

# GETRF
	# compare to dense solve

	

import numpy as np

import hMatrix
import hMatrix_todense
import hMatrixGEMM
import hMatrixTRSM
import hMatrixGETRF


def rand_complex( shp ):
	return 2*np.random.random( shp )-1 + 1j*( 2*np.random.random( shp )-1 )


def make_random_hMatrix( M, N, block_size, dense_diag=True ):
	
	tru_A = np.zeros( (M,N), dtype=np.complex128 )
	A = hMatrix.hMatrix( M, N, block_size )
	
	job_stack = [( A.root_node(), 0,M, 0,N )]
	while len(job_stack) > 0 :

		node, grb,gre, gcb,gce = job_stack.pop(0)
		m,n = node.shape
		on_diag = (grb,gre) == (gcb,gce)
		n_row_blocks = hMatrix.num_blocks( block_size, m )
		n_col_blocks = hMatrix.num_blocks( block_size, n )
		
		if on_diag and dense_diag:
			if n_row_blocks == 1 :
				# Make dense
				D = rand_complex( (m,n) )
				tru_A[ grb:gre, gcb:gce ] = D[:,:]
				node.insert_dense( D )
			
			else: # on diag and can partition
				# Dense or partition
				dense_prob = 0.1 if n_row_blocks > 3 else 0.3
				
				r = np.random.rand()
				
				if r < dense_prob :
					# Make dense
					D = rand_complex( (m,n) ) + 2*np.eye( m, dtype=np.complex128 )
					tru_A[ grb:gre, gcb:gce ] = D[:,:]
					node.insert_dense( D )
				
				else:
					# Partition
					bounds = [ hMatrix.block_begin( block_size, m, i ) for i in range( 1, n_row_blocks, 1 ) ]
					n_part_pts = np.random.randint( 1, min( len(bounds), 4 )+1 )
					part_pts = np.random.choice( bounds, size=n_part_pts, replace=False )
					part_pts.sort()
					
					node.partition( part_pts, part_pts )
					bounds = node.get_row_bounds()
					
					for i in range(node.children_shape()[0]):
						for j in range(node.children_shape()[1]):
							child = node.get_child( i,j )
							rb = grb + bounds[i]
							re = grb + bounds[i+1]
							cb = gcb + bounds[j]
							ce = gcb + bounds[j+1]
							job_stack.append(( child, rb,re, cb,ce ))
					
		
		else :
			if n_row_blocks == 1 and n_col_blocks == 1 :
				
				r = np.random.rand()
				
				if r < 0.25 :
					# Make dense
					D = rand_complex( (m,n) )
					tru_A[ grb:gre, gcb:gce ] = D[:,:]
					node.insert_dense( D )
				
				elif r < 0.80 :
					# Make lowrank
					rank = np.random.randint( 3,6 )
					L = rand_complex( (m,rank) )
					R = rand_complex( (rank,n) )
					tru_A[ grb:gre, gcb:gce ] = L @ R
					node.insert_lowrank( L, R )
			
			else : # off diag and can partition
				# dense or zero or lowrank or partition

				r = np.random.rand()
				
				if r < 0.05 :
					# Make dense
					D = rand_complex( (m,n) )
					tru_A[ grb:gre, gcb:gce ] = D[:,:]
					node.insert_dense( D )
				
				elif r < 0.50 :
					# Make lowrank
					rank = np.random.randint( 3,6 )
					L = rand_complex( (m,rank) )
					R = rand_complex( (rank,n) )
					tru_A[ grb:gre, gcb:gce ] = L @ R
					node.insert_lowrank( L, R )
				
				elif r < 0.95:
					# Partition
					if n_row_blocks > 1 :
						row_bounds = [ hMatrix.block_begin( block_size, m, i ) for i in range( 1, n_row_blocks, 1 ) ]
						n_row_part_pts = np.random.randint( 1, min( len(row_bounds), 4 )+1 )
						row_part_pts = np.random.choice( row_bounds, size=n_row_part_pts, replace=False )
						row_part_pts.sort()
					else:
						row_part_pts=[]
					
					if n_col_blocks > 1 :
						col_bounds = [ hMatrix.block_begin( block_size, n, i ) for i in range( 1, n_col_blocks, 1 ) ]
						n_col_part_pts = np.random.randint( 1, min( len(col_bounds), 4 )+1 )
						col_part_pts = np.random.choice( col_bounds, size=n_col_part_pts, replace=False )
						col_part_pts.sort()
					else:
						col_part_pts=[]
					
					node.partition( row_part_pts, col_part_pts )
					row_bounds = node.get_row_bounds()
					col_bounds = node.get_col_bounds()
					
					for i in range(node.children_shape()[0]):
						for j in range(node.children_shape()[1]):
							child = node.get_child( i,j )
							rb = grb + row_bounds[i]
							re = grb + row_bounds[i+1]
							cb = gcb + col_bounds[j]
							ce = gcb + col_bounds[j+1]
							job_stack.append(( child, rb,re, cb,ce ))

	return A,tru_A


def test_construction():

	seed = np.random.randint(0,10000)
#	seed = 2938
	print(f"Random seed is {seed}")
	np.random.seed(seed)

	A, tru_dense_A = make_random_hMatrix( 43,43, 5 )
	
	dense_A = hMatrix_todense.hMatrix_todense( A )
	
	assert np.max( np.abs( tru_dense_A - dense_A ) ) < 1e-9
	
	n_blocks = hMatrix.num_blocks( 5, 43 )
	bounds = [ hMatrix.block_begin( 5, 43, i ) for i in range(n_blocks+1) ]
	for i in range(n_blocks):
		for j in range(n_blocks):
			rb = bounds[i]
			re = bounds[i+1]
			cb = bounds[j]
			ce = bounds[j+1]
			
			tmp = hMatrix_todense.hMatrix_todense( A[rb:re,cb:ce] )
			assert np.max( np.abs( tru_dense_A[rb:re,cb:ce] - tmp ) ) < 1e-9


	




def test_hMatrixGEMM_1():

	seed = np.random.randint(0,10000)
	print(f"Random seed is {seed}")
	np.random.seed(seed)

	A = hMatrix.hMatrix( 20,14, 5 )
	B = hMatrix.hMatrix( 14, 4, 5 )
	C = hMatrix.hMatrix( 20, 4, 5 )
	root = A.root_node()
	root.partition( [5,10], [6] )
	cshape = root.children_shape()
	for i in range(cshape[0]):
		for j in range(cshape[1]):
		    child = root.get_child(i,j)
		    shp = child.shape
		    child.insert_dense( np.astype(np.random.random(shp),np.complex128) )
	a = hMatrix_todense.hMatrix_todense( A )

	root = B.root_node()
	shp = root.shape
	rank = 3
	root.insert_lowrank(  np.astype(np.random.random((shp[0],rank)),np.complex128),  np.astype(np.random.random((rank,shp[1])),np.complex128) )
	b = hMatrix_todense.hMatrix_todense( B )
	hMatrixGEMM.hMatrixGEMM( 1.0, A, B, C )
	c = hMatrix_todense.hMatrix_todense( C )

	assert np.linalg.norm((a@b)-c) < 1e-9


#def test_hMatrixGEMM_2(): # test support for size 0 matrices
#	
#	A = hMatrix.hMatrix( 10,10,10 )
#	A.insert_dense( A.root_node(), np.eye( 10 , dtype=A.dtype ) )
#	B = hMatrix.hMatrix( 0,1,10 )
#	C = hMatrix.hMatrix( 0,1,10 )
#	
#	hMatrixGEMM.hMatrixGEMM( 1.0, A[1:1,1:1], B, C )
	

def test_hMatrixGEMM_2():

	seed = np.random.randint(0,10000)
	print(f"Random seed is {seed}")
	np.random.seed(seed)

	A, tru_dense_A = make_random_hMatrix( 43,63, 5 )
	B, tru_dense_B = make_random_hMatrix( 63,43, 5 )
	C, tru_dense_C = make_random_hMatrix( 43,43, 5 )
	
	hMatrixGEMM.hMatrixGEMM( -1.0, A, B, C )
	
	dense_C = hMatrix_todense.hMatrix_todense( C )
	
	tru_dense_C -= tru_dense_A @ tru_dense_B
	
	assert np.max( np.abs( tru_dense_C - dense_C ) ) < 1e-9



from scipy.linalg import solve_triangular

def test_hMatrixTRSM_1():

	seed = np.random.randint(0,10000)
	seed = 6063
	print(f"Random seed is {seed}")
	np.random.seed(seed)
	
	M = 43
	block_size = 5

	A, tru_dense_A = make_random_hMatrix( M, M, block_size )
	
	p = np.zeros( M, dtype=np.uint64 )
	n_blocks = hMatrix.num_blocks( block_size, M )
	for i in range(n_blocks):
		begin = hMatrix.block_begin( block_size, M, i )
		end = hMatrix.block_begin( block_size, M, i+1 )
		p[begin:end] = np.arange(end-begin)
	
	
	
	B, tru_dense_B = make_random_hMatrix( M, 10, block_size )
	
	side = "L"
	uplo = "L"
	hMatrixTRSM.hMatrixTRSM( side, uplo, p, A, B )
	
	dense_x = hMatrix_todense.hMatrix_todense( B )
	
	tru_dense_x = solve_triangular( tru_dense_A, tru_dense_B, trans=0, lower=True, unit_diagonal=True, overwrite_b=False, check_finite=False)
	
	assert np.max( np.abs( tru_dense_x - dense_x ) ) < 1e-9
	
	
	
	B, tru_dense_B = make_random_hMatrix( 10, M, block_size )
	
	side = "R"
	uplo = "L"
	hMatrixTRSM.hMatrixTRSM( side, uplo, p, A, B )
	
	dense_x = hMatrix_todense.hMatrix_todense( B )
	
	tru_dense_x = solve_triangular( tru_dense_A, tru_dense_B.T, trans=1, lower=True, unit_diagonal=True, overwrite_b=False, check_finite=False)
	
	assert np.max( np.abs( tru_dense_x.T - dense_x ) ) < 1e-9
	
	
	
	B, tru_dense_B = make_random_hMatrix( M, 10, block_size )
	
	side = "L"
	uplo = "U"
	hMatrixTRSM.hMatrixTRSM( side, uplo, p, A, B )
	
	dense_x = hMatrix_todense.hMatrix_todense( B )
	
	tru_dense_x = solve_triangular( tru_dense_A, tru_dense_B, trans=0, lower=False, unit_diagonal=False, overwrite_b=False, check_finite=False)
	
	assert np.max( np.abs( tru_dense_x - dense_x ) ) < 1e-6
	

	
	B, tru_dense_B = make_random_hMatrix( 10, M, block_size )
	
	side = "R"
	uplo = "U"
	hMatrixTRSM.hMatrixTRSM( side, uplo, p, A, B )
	
	dense_x = hMatrix_todense.hMatrix_todense( B )
	
	tru_dense_x = solve_triangular( tru_dense_A, tru_dense_B.T, trans=1, lower=False, unit_diagonal=False, overwrite_b=False, check_finite=False)
	
	assert np.max( np.abs( tru_dense_x.T - dense_x ) ) < 1e-6


















