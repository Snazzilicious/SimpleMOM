
import numpy as np

def Leaf_todense( A, D ):
	if A.block_type() == "Dense" :
		D_a = A.get_dense_data()
		D[:,:] = D_a[:,:]
	elif A.block_type() == "LowRank" :
		L_a, R_a = A.get_lowrank_data()
		D[:,:] = L_a @ R_a


def queue_H_todense( A, D ):
	row_bounds = A.get_row_bounds()
	col_bounds = A.get_col_bounds()
	
	jobs=[]
	for row_begin,row_end in zip(row_bounds[:-1],row_bounds[1:]):
		for col_begin,col_end in zip(col_bounds[:-1],col_bounds[1:]):
			a = A[ row_begin:row_end, col_begin:col_end ]
			d = D[ row_begin:row_end, col_begin:col_end ]
			
			jobs.append(( a,d ))
	
	return jobs


def hMatrix_todense( A ):
	D = np.zeros( A[:,:].shape, dtype=A.dtype )
	
	job_stack = [( A[:,:], D[:,:] )]
	while len(job_stack) > 0 :

		a,d = job_stack.pop(0)
		if a.shape != d.shape:
			raise AssertionError(f"Incompatible matrix dimensions: {a.shape}, {d.shape}")
		
		type_a = a.block_type()
		
		# Trivial case
		if type_a == "Zero" :
			continue
		
		# Base case
		elif type_a != "H" :
			Leaf_todense( a, d )

		# General case
		else :
			job_stack[0:0] = queue_H_todense( a, d )
	
	return D
