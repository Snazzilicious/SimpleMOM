

import hMatrix
import hMatrix_todense

import numpy as np

# TODO
# check node types and dimensions and bounds
# try getting matrices
# compare to todense
# check slice types and dimensions and bounds
# try getting matrices
# compare to todense

# GEMM
	# rSVD x2
	# llaxpy
	# leaf gemm
	# GEMM - compare todense multiplies

# TRSM
	# compare to todense tri solve

# GETRF
	# compare to dense solve


def test_construction():
	
	A = hMatrix.hMatrix( 20,20, 5 )
	root = A.root_node()
	A.partition( root, [0,10,20], [0,10,20] )
	
	assert root.shape == (20,20)
	assert root.block_type() == "H"
	
	for i in range(2):
		for j in range(2):
		    if not (i == 0 and j == 0) :
		        root.get_child(i,j).insert_dense( np.eye(10,dtype=np.complex128) )

	child = root.get_child( 0,0 )
	child.partition([0,5],[5,10])

	for i in range(2):
		for j in range(2):
		    child.get_child(i,j).insert_dense( np.eye(5,dtype=np.complex128) )

	x = hMatrix_todense.hMatrix_todense( A )



