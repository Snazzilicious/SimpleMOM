
import hMatrix

import numpy as np

def test_num_blocks():
	
	assert hMatrix.num_blocks( 10, 23 ) == 2
	
	assert hMatrix.num_blocks( 10, 5 ) == 1
	
	assert hMatrix.num_blocks( 256, 2e6 ) == 7812


def test_block_size():
	
	min_block_size = 200
	length = 34567
	all_sizes = [hMatrix.block_size( min_block_size, length, i ) for i in range(hMatrix.num_blocks( min_block_size, length ))]
	
	assert length == sum(all_sizes)
	
	min_block_size = 200
	length = 100
	all_sizes = [hMatrix.block_size( min_block_size, length, i ) for i in range(hMatrix.num_blocks( min_block_size, length ))]
	
	assert length == sum(all_sizes)


def test_block_begins():
	
	min_block_size = 200
	length = 368897
	all_sizes = np.array([hMatrix.block_size( min_block_size, length, i ) for i in range(hMatrix.num_blocks( min_block_size, length ))])
	
	bounds = [ hMatrix.block_begin( min_block_size, length, i ) for i in range(hMatrix.num_blocks( min_block_size, length )+1) ]
	
	assert bounds[0] == 0
	assert bounds[-1] == length
	assert np.all( all_sizes == np.diff(bounds) )


def test_block_index():
	
	min_block_size = 20
	length = 369
	bounds = [ hMatrix.block_begin( min_block_size, length, i ) for i in range(hMatrix.num_blocks( min_block_size, length )+1) ]
	
	for i in range(len(bounds)-1):
		begin = bounds[i]
		end = bounds[i+1]
		
		inds = np.array([ hMatrix.block_index( min_block_size, length, i ) for i in range(begin,end,1) ])
		
		assert np.all( inds == i )
		

def test_on_block_boundary():
	
	min_block_size = 20
	length = 369
	bounds = [ hMatrix.block_begin( min_block_size, length, i ) for i in range(hMatrix.num_blocks( min_block_size, length )+1) ]
	
	for i in range(length+1):
		assert hMatrix.on_block_boundary( min_block_size, length, i ) == (i in bounds)
	
	
	
