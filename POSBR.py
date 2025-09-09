

import numpy as np
import ctypes

nanort_lib = ctypes.CDLL( "RayTracing/nanortlib.so" )
nanort_lib.get_handle.restype = ctypes.c_void_p
nanort_lib.intersects_first_interface.restype = None
nanort_lib.free_handle.restype = None


class RayTracer:
	def __init__( self, vertices, faces ):
		
		n_faces = len(faces)
		vertices_ptr = vertices.ctypes.data_as( ctypes.POINTER(ctypes.c_double) )
		faces_ptr = faces.ctypes.data_as( ctypes.POINTER(ctypes.c_uint) )
		
		self.nanort_handle = nanort_lib.get_handle( n_faces, vertices_ptr, faces_ptr )
	
	def intersects_first( self, ray_orgs, ray_dirs ):
		
		n_rays = len(ray_orgs)
		hit_indices = -np.ones( n_rays, dtype=np.int64 ) # consistent with c_long according to np.ctypeslib.as_ctypes_type(np.int64)
		
		orgs_ptr = ray_orgs.ctypes.data_as( ctypes.POINTER(ctypes.c_double) )
		dirs_ptr = ray_dirs.ctypes.data_as( ctypes.POINTER(ctypes.c_double) )
		hits_ptr = hit_indices.ctypes.data_as( ctypes.POINTER(ctypes.c_long) )
		
		nanort_lib.intersects_first_interface( n_rays, orgs_ptr, dirs_ptr, hits_ptr, self.nanort_handle )
		
		return hit_indices
	
	
	def __del__(self):
		nanort_lib.free_handle( self.nanort_handle )
