
import os
import sys

parent_dir = os.path.sep.join( os.path.abspath(__file__).split(os.path.sep)[:-2] )
sys.path.insert( 0, parent_dir )

import Geometry as geo

import numpy as np



def is_unit_vector( v ):
	return np.isclose( np.linalg.norm(v), 1.0 )

def test_to_cartesian():
	
	# Easy ones
	rs = [ 0.0, 0.6, 1.1 ]
	x_angles = [ i*np.pi/6 for i in [1, 3, 4, 7, 10] ]
	
	r = 1.1
	x_angle = np.pi / 6
	z_angle = 0.0
	
	# r is zero
	
	# z_angle is zero or \pi
	
	# mismatching dimensions in inputs

	assert True


def test_to_spherical():

	# Easy ones
	
	# x and y are zero
	
	# y and z are zero
	
	# x, y, and z are zero
	
	# mismatching dimensions in inputs

	assert True

def test_unit_vectors():
	
	# Easy ones
	
	# z_angle is zero or \pi
	
	# variety of dimensions of inputs
	
	assert True
