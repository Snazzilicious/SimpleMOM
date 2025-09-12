"""Functions for reading and writing outputs from simulations.

Conventions:
	Azimuthal angle is referred to as 'X angle'
	Zenith angle is referred to as 'Z angle'
	Scattered fields transmit and receive angles are given in the order transmit x, transmit z, receive x, receive z
	Angles are in radians
"""
import numpy as np

from scipy.constants import c

def aurora_field_scale_factor( freq ):
	"""Aurora normalizes incident wave to 1 V / m, so fields need to be scaled by 1/\lambda
	
	Arguments
		freq : float
			The frequency of the incident wave in Hertz.
	
	Returns
		aurora_field_scale_factor : float
			Factor by which to multiply Aurora's output fields
			to become consistent with a unit polarization incident field.
	"""
	return freq / c

import h5py

def get_aurora_field( filename ):
	"""Retrieves Aurora output fields from its output h5 file.
	
	Arguments
		filename : string
			Name of the .h5 file from which to get scattered fields.
	
	Returns
		results : dictionary
			Fields from the .fld file organized by polarization in and out. Keys are "ZZ", "ZX", "XZ", or "XX".
			Values are tuples where the first entry is an ndarray of shape ( num_in_and_out_angles, 4 ) 
			of the the incident and scattered angles in radians, and the second entry is a numpy array
			of length num_in_and_out_angles of the complex scattered field components.
	"""
	f = h5py.File(filename, 'r')
	hash_string = list(f['aurorae']['results'].keys())[0]
	
	xz_in_angs = np.vstack([[j for j in i] for i in list(f['aurorae']['results'][hash_string]['excitations']['planewave']['table'])])
	xz_in_angs = xz_in_angs[:,1:]
	
	dat = np.vstack([[j for j in i] for i in list(f['aurorae']['results'][hash_string]['observations']['scattered farfield']['table'])])
	
	xz_out_angles = dat[:,1:3]
	xzxz_angles = np.column_stack(( xz_in_angs, xz_out_angles ))
	
	freq = dat[:,3]
	
	results = {}
	results["ZZ"] = ( xzxz_angles, dat[:,4]+1j*dat[:,5] )
	results["ZX"] = ( xzxz_angles, dat[:,8]+1j*dat[:,9] )
	results["XZ"] = ( xzxz_angles, dat[:,10]+1j*dat[:,11] )
	results["XX"] = ( xzxz_angles, dat[:,6]+1j*dat[:,7] )
	# TODO divide by scale
	
	return results # freq


def stars_field_scale_factor():
	"""Stars normalizes incident wave so that H0 = 1 A/m, or E0 = 376 V/m, so fields need to be scaled by 1/376.
	
	Returns
		stars_field_scale_factor : float
			Factor by which to multiply Stars's output fields
			to become consistent with a unit polarization incident field.
	"""
	return 1.0/376.0

import regex as re

def get_stars_fld( filename ):
	"""Retrieves Stars output fields from its .fld file.
	
	Arguments
		filename : string
			Name of the .fld file from which to get scattered fields.
	
	Returns
		results : dictionary
			Fields from the .fld file organized by polarization in and out. Keys are "ZZ", "ZX", "XZ", or "XX".
			Values are tuples where the first entry is an ndarray of shape ( num_in_and_out_angles, 4 ) 
			of the the incident and scattered angles in radians, and the second entry is a numpy array
			of length num_in_and_out_angles of the complex scattered field components.
	"""
	f=open(filename , "r" )
	for _ in range(10):
		f.readline()
	
	line = f.readline()
	m = re.findall( "Re\[([TP]{2})\]", line )
	key_change = { "TT" : "ZZ", "TP" : "ZX", "PT" : "XZ", "PP" : "XX" }
	keys = [ key_change[i] for i in m ]
	
	f.readline()
	
	dat=[]
	for line in f:
		dat.append([float(v) for v in line.split(",")])
	f.close()
	dat = np.vstack(dat)

	xzxz_angs = np.column_stack(( dat[:,1], dat[:,0], dat[:,3], dat[:,2] ))
	xzxz_angs = np.deg2rad( xzxz_angs )
	
	results = { k : ( xzxz_angs, dat[:,4+2*i] + 1j*dat[:,4+2*i+1] ) for i,k in enumerate(keys) }

	return results

# TODO
def get_stars_current( filename ):
	f = open(filename,"r")
	f.readline()
	
	coord=[]
	dat1=[]
	dat2=[]
	for line in f:
		tmp = [ float(s) for s in line.split()[1:] ]
		coord.append(tmp[:3])
		dat1.append([ complex(tmp[3],tmp[4]), complex(tmp[5],tmp[6]), complex(tmp[7],tmp[8]) ])
		dat2.append([ complex(tmp[9],tmp[10]), complex(tmp[11],tmp[12]), complex(tmp[13],tmp[14]) ])
	coords.append(np.vstack(coord))
	Js1.append(np.vstack(dat1))
	Js2.append(np.vstack(dat2))
	
	f.close()



import regex as re

def gems_field_scale_factor():
	return 1.0/3.545

def get_gems_field( filename ):
	"""Retrieves GEMS output fields from its .rcs file.
	
	Arguments
		filename : string
			Name of the .rcs file from which to get scattered fields.
	
	Returns
		results : dictionary
			Fields from the .fld file organized by polarization in and out. Keys are "ZZ", "ZX", "XZ", or "XX".
			Values are tuples where the first entry is an ndarray of shape ( num_in_and_out_angles, 4 ) 
			of the the incident and scattered angles in radians, and the second entry is a numpy array
			of length num_in_and_out_angles of the complex scattered field components.
	"""
	f = open(filename,"r")
	
	key_change = { "VV" : "ZZ", "VH" : "ZX", "HV" : "XZ", "HH" : "XX" }

	run_id = 0
	runs = []
	
	mode = "FIND_RUN"
	for line in f:
		if mode == "FIND_RUN":
			m = re.search( "Run # *([0-9]+)", line )
			if m is not None:
				if int(m[1]) != run_id:
					run_id = int(m[1])
					runs.append({})
				mode = "FIND_START"
		
		elif mode == "FIND_START":
			if 'Start RCS Data' in line:
				mode = "COLLECT"
		
		elif mode == "COLLECT":
			results = runs[-1]
			
			if 'End RCS Data' in line:
				mode = "FIND_RUN"
				
			else:
				spline = line.split()
				key = key_change[ spline[2][1]+spline[5][1] ]
				xzxz_angs = [ float(i) for i in [ spline[0], spline[1], spline[3], spline[4] ] ]
				fld = float(spline[6])+1j*float(spline[7])
				
				if not key in results.keys():
					results[key] = ([],[])
			
				results[ key ][0].append( xzxz_angs )
				results[ key ][1].append( fld )
				
	f.close()
	
	for results in runs:
		for k in results:
			adj_xzxz = np.deg2rad(results[k][0])
			adj_xzxz[:,1] *= -1.0
			adj_xzxz[:,1] += np.pi/2
			adj_xzxz[:,3] *= -1.0
			adj_xzxz[:,3] += np.pi/2
			
			results[k] = ( adj_xzxz, np.array(results[k][1]) )
	
	
	if len(runs) == 1:
		return runs[0]
	else:
		return runs


from sklearn.neighbors import KNeighborsRegressor
from EM_Utils import x_angle_z_angle_to_dir

class FieldRegressor :
	"""Function like object used to evaluate output fields at various incident and scattered angles.
	"""
	def xzxz2xyzxyz( xzxz_angs ):
		"""Get cartesian coordinates on the unit sphere from x & z angles.
		The regressor is based on Euclidean distance and so spherical coordinates are problematic.
		"""
		xyzxyz = np.zeros((len(xzxz_angs),6))
		xyzxyz[:,:3] = x_angle_z_angle_to_dir( xzxz_angs[:,0], xzxz_angs[:,1] )
		xyzxyz[:,3:] = x_angle_z_angle_to_dir( xzxz_angs[:,2], xzxz_angs[:,3] )
		
		return xyzxyz
	
	def __init__( self, xzxz_angs, fld ):
		"""Construct and fit the object to field data.
		
		Arguments
			xzxz_angs : ndarray of shape ( num_samples, 4 )
				Incident and scattered angles for each field sample.
			fld : numpy array of complex of shape ( num_samples[, ndim] )
				Scattered field to which to fit regressor.
		"""
		xyzxyz = FieldRegressor.xzxz2xyzxyz( xzxz_angs )
		
		self.feature_shape = fld.shape[1:]
		self.complex_field = np.any(np.iscomplex(fld))
		
		if self.complex_field:
			fld = np.column_stack(( np.real(fld), np.imag(fld) ))
		
		self.knr = KNeighborsRegressor( n_neighbors=1, weights='distance' )
		self.knr.fit( xyzxyz, fld )
	
	def __call__( self, xzxz_angs ):
		"""Evaluate scattered field at specified points.
		
		Arguments
			xzxz_angs : ndarray of shape ( num_eval_pts, 4 )
				Incident and scattered angles at which to calculate the field.
			
		Returns
			fld : numpy array of complex of shape ( num_eval_pts[, ndim] )
				Scattered field computed at each evaluation point.
		"""
		xyzxyz = FieldRegressor.xzxz2xyzxyz( xzxz_angs )
		
		pred = self.knr.predict( xyzxyz )
		
		if self.complex_field:
			ncol = pred.shape[1]
			pred = pred[:,:int(ncol/2)] + 1j*pred[:,int(ncol/2):]
		
		return pred.reshape( (-1,) + self.feature_shape )
		


