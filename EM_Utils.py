"""Formulas from the field of Electromagnetic modeling.
"""


import numpy as np

from scipy.constants import c, epsilon_0, mu_0
eps = epsilon_0
mu = mu_0

def omega_from_Hz( freq ):
	"""Compute the angular frequency in the time-dependent part of time-harmonic solutions e^{i \omega t}.
	
	Arguments
		freq : numpy array of floats
			The frequencies in Hertz.
	
	Returns
		omega : numpy array of floats
			The angular frequencies in radians per second.
	"""
	return 2*np.pi*freq

def k_from_omega( w ):
	"""Compute the wave number.
	
	Arguments
		omega : numpy array of floats
			The angular frequencies in radians per second.
	
	Returns
		k : numpy array of floats
			The wave numbers in radians per meter.
	"""
	return w/c

def omega_from_k( k ):
	"""Compute the angular frequency in the time-dependent part of time-harmonic solutions e^{i \omega t}.
	
	Arguments
		k : numpy array of floats
			The wave numbers in radians per meter.
	
	Returns
		omega : numpy array of floats
			The angular frequencies in radians per second.
	"""
	return c*k


"""Plane wave E-H conversions
"""

def H_from_E( p,E ):
	"""Computes the H vector associated with a plane wave with given direction and E vector.
	
	Arguments
		p : ...
			The propagation direction of the plane wave.
		E : ...
			The electric field polarization vector of the plane wave.
	
	Returns
		H : ...
			The magnetic polarization vector of the plane wave.
	"""
	return np.sqrt( eps/mu ) * np.cross( p, E )

def E_from_H( p,H ):
	"""Computes the E vector associated with a plane wave with given direction and H vector.
	
	Arguments
		p : ...
			The propagation direction of the plane wave.
		H : ...
			The magnetic polarization vector of the plane wave.
	
	Returns
		E : ...
			The electric field polarization vector of the plane wave.
	"""
	return np.sqrt( mu/eps ) * np.cross( H, p )


def planewave( k, prop, pol, pts ):
	"""Functional form of a plane wave.
	
	Arguments
	
	Returns
	"""
	return np.outer( np.exp( -1j * k * pts.dot(prop) ), pol )


"""Green's function
"""

def G( x, y, k ):
	"""Helmholtz Green's function.
	
	Arguments
		x : array-like of floats
			The field or observation point.
		y : array-like of floats
			The source point.
		k : float
			The wave number in radians per meter.
	
	Returns
		G : array-like of floats
			The value of the Helmholtz Green's function evaluted at all pairs of points passed in.
	"""
	R = np.linalg.norm( x-y, axis=1 )
	return -np.exp( -1j*k*R ) / ( 4*np.pi * np.maximum( R, 1e-15 ) )

def gradx_G( x, y, k ):
	"""Gradient of the Helmholtz Green's function with respect to x.
	
	Arguments
		x : array-like of floats
			The field or observation point.
		y : array-like of floats
			The source point.
		k : float
			The wave number in radians per meter.
	
	Returns
		gradG : array-like of floats
			The gradient of the Helmholtz Green's function w.r.t. x, evaluted at all pairs of points passed in.
	"""
	r = x - y
	R = np.linalg.norm( r, axis=1 )
	return ( r.T * ( 1j*k*R + 1 ) * np.exp( -1j*k*R ) / ( 4*np.pi * np.maximum( R**3, 1e-15 ) ) ).T


def G_ff( p, y, k ):
	# TODO use this in farfield functions - need v&v tests before updating
	return np.exp( 1j * k * np.sum( pts * p, axis=1 ) ) / (4*np.pi)




"""Working with E field
"""
def PO_current_from_E( prop_dir, E_pol, normals ):
	return 2 * np.sqrt( eps/mu ) * np.cross( normals, np.cross( prop_dir, E_pol ) )

def IPO_current_from_E( prop_dir, E_pol, normals, tol=1e-1 ):
	illuminated = np.sum( normals * prop_dir, axis=1 ) < -tol
	currents = PO_current_from_E( prop_dir, E_pol, normals )
	return ( currents.T * illuminated ).T

def E_farfield( obs_dir, k, currents, pts, wts ):
	w = omega_from_k(k)
	
	ff = ( 1j * w * mu / (4*np.pi) ) * np.sum( wts * np.exp( 1j*k*pts.dot(obs_dir) ) * currents.T, axis=1 )
	ff -= ff.dot(obs_dir) * obs_dir
	
	return ff

def bistatic_E_field( x_z_angles, k, currents, pts, wts ):
	ffs=np.zeros([len(x_z_angles),3],dtype=np.complex128)
	
	for i,(x_angle,z_angle) in enumerate(x_z_angles):
		obs_dir = x_angle_z_angle_to_dir( x_angle, z_angle )
		ffs[i] = E_farfield( obs_dir, k, currents, pts, wts )
	
	return ffs

def E_farfield_at_pos( origin, pos, k, currents, pts, wts ):
	ffs=np.zeros([len(pos),3],dtype=np.complex128)
	
	for i,(o,p) in enumerate(zip(origin,pos)):
		R = p-o
		norm_R = np.linalg.norm(R)
		obs_dir = R/norm_R
		ffs[i] = np.exp( -1j * k * norm_R ) * E_farfield( obs_dir, k, currents, pts-o, wts ) / norm_R
	
	return ffs


"""Working with H field
"""
def PO_current_from_H( H_pol, normals ):
	return 2 * np.cross( normals, H_pol )

def IPO_current_from_H( prop_dir, H_pol, normals, tol=1e-1 ):
	illuminated = np.sum( normals * prop_dir, axis=1 ) < -1e-1
	currents = PO_current_from_H( H_pol, normals )
	return ( currents.T * illuminated ).T

def bistatic_H_field( x_z_angles, k, currents, pts, wts ):
	
	ffs=np.zeros([len(x_z_angles),3],dtype=np.complex128)
	
	for i,(x_angle,z_angle) in enumerate(x_z_angles):
		obs_dir = x_angle_z_angle_to_dir( x_angle, z_angle )
		V = np.cross( currents, obs_dir )
		ffs[i] = ( -1j * k / (4*np.pi) ) * np.sum( wts * np.exp( 1j*k*pts.dot(obs_dir) ) * V.T, axis=1 )
	
	return ffs

def bistatic_H_field_contribution( x_z_angles, k, currents, pts, wts ):
	
	contribs=[]
	
	for i,(x_angle,z_angle) in enumerate(x_z_angles):
		obs_dir = x_angle_z_angle_to_dir( x_angle, z_angle )
		V = np.cross( currents, obs_dir )
		contribs.append( ( -1j * k / (4*np.pi) ) * wts * np.exp( 1j*k*pts.dot(obs_dir) ) * V.T )
	
	return contribs


""" For presenting farfield """

def projected_farfield( x_z_angles, ffs ):
	"""Project farfield onto spherical surface basis vectors.
	
	Arguments
		x_z_angles : ndarray of float of shape ( num_field_samples, 2 )
			The scattering angles of the farfield samples.
		ffs : ndarray of complex of shape ( num_field_samples, 3 )
			The farfield vector in cartesian 3 space to be projected onto
			the spherical surface basis vectors.
	
	Returns
		coefs_x_angle : numpy array of complex of length num_field_samples
			The farfield component in the x-angle direction.
		coefs_z_angle : numpy array of complex of length num_field_samples
			The farfield component in the z-angle direction.
	"""
	ref_vecs_z = z_angle_hat( x_z_angles[:,0], x_z_angles[:,1] )
	ref_vecs_x = x_angle_hat( x_z_angles[:,0], x_z_angles[:,1] )
	coefs_z_angle = np.sum( ffs * ref_vecs_z, axis=1 )
	coefs_x_angle = np.sum( ffs * ref_vecs_x, axis=1 )

	return coefs_x_angle, coefs_z_angle


def to_dBsm( x ):
	"""Convert to decibels per square meter.
	For RCSs, pass in |E_s|/|E_i|.
	
	Arguments
		x : ndarray of complex
			Value(s) to be converted.
	
	Returns
		dBsm : ndarray of floats
			10*log_10( 4\pi (x)^2 + 1e-12 )
	"""
	# For RCS: 10*log_10( 4\pi ( |E_s|/|E_i| )^2 )
	return 10*np.log10( 4*np.pi*np.abs( x )**2+1e-12 )





