"""Geometric mathematical functions and formulas.
"""

import numpy as np


"""Spherical coordinates functions
"""

def cartesian_to_spherical( x, y, z ):
	"""Convert cartesian coordinates to spherical coordinates.
	
	Arguments
		x : # XXX look at doc strings for called functions XXX
			Distance along the x-axis.
		y : 
			Distance along the y-axis.
		z : 
			Distance along the z-axis.
	
	Returns
		r : 
			Distance from the origin.
		x_angle : 
			Angle from the x-axis in the XY plane, measured in radians.
		z_angle : 
			Angle from the z-axis, measured in radians.
	"""
	r = np.sqrt( x**2 + y**2 + z**2 )
	
	x_angle = np.arctan2( y, x )
	
	z_angle = np.arccos( z / (r+1e-12) )
	z_angle[r < 1e-12] = 0.0
	
	return r, x_angle, z_angle


def spherical_to_cartesian( r, x_angle, z_angle ):
	"""Convert spherical coordinates to cartesian coordinates.
	
	Arguments
		r : 
			Distance from the origin.
		x_angle : 
			Angle from the x-axis in the XY plane, measured in radians.
		z_angle : 
			Angle from the z-axis, measured in radians.
	
	Returns
		x : 
			Distance along the x-axis.
		y : 
			Distance along the y-axis.
		z : 
			Distance along the z-axis.
	"""
	
	x = r * np.cos(x_angle) * np.sin(z_angle)
	y = r * np.sin(x_angle) * np.sin(z_angle)
	
	rx_shape = np.broadcast_shapes( np.shape(r), np.shape(x_angle) )
	z = np.broadcast_to( r, rx_shape ) * np.cos(z_angle)
	
	return x, y, z


def r_hat( x_angle, z_angle ):
	"""Unit vector in the spherical coordinate 'r' direction.
	
	Arguments
		x_angle : 
			Angle from the x-axis in the XY plane, measured in radians.
		z_angle : 
			Angle from the z-axis, measured in radians.
	
	Returns
		r_hat : 
			The unit vector in the coordinate 'r' direction at the given angles.
	"""
	dirs = np.column_stack( spherical_to_cartesian( 1.0, x_angle, z_angle ) )
	return np.squeeze( dirs )

def x_angle_hat( x_angle, z_angle=None ):
	"""Unit vector in the spherical coordinate 'x_angle' direction.
	
	Arguments
		x_angle : 
			Angle from the x-axis in the XY plane, measured in radians.
		z_angle : 
			Angle from the z-axis, measured in radians.
	
	Returns
		x_angle_hat : 
			The unit vector in the coordinate 'x_angle' direction at the given angles.
	"""
	return r_hat( x_angle+np.pi/2, np.pi/2 )

def z_angle_hat( x_angle, z_angle ):
	"""Unit vector in the spherical coordinate 'z_angle' direction.
	
	Arguments
		x_angle : 
			Angle from the x-axis in the XY plane, measured in radians.
		z_angle : 
			Angle from the z-axis, measured in radians.
	
	Returns
		z_angle_hat : 
			The unit vector in the coordinate 'z_angle' direction at the given angles.
	"""
	return -r_hat( x_angle+np.pi, np.pi/2 - z_angle )

def x_angle_z_angle_to_dir( x_angle, z_angle ):
	"""Alias for `r_hat`, converts direction given by spherical coordinates angles to
	a cartesian unit vector.
	"""
	return r_hat( x_angle, z_angle )
