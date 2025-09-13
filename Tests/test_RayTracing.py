
import os
import sys

parent_dir = os.path.sep.join( os.path.abspath(__file__).split(os.path.sep)[:-2] )
sys.path.insert( 0, parent_dir )

import POSBR
import numpy as np

def get_test_data():
	vertices = np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0]],dtype=np.float64)
	faces = np.array([[0,1,2],[1,3,2]],dtype=np.int32)

	orgs = np.array([[0.25,0.25,0.33],[0.4,0.5,1],[0.5,0.75,1.5],[0.5,0.75,1],[0.2,0.2,0.2]],dtype=np.float64)
	dirs = -np.array([[0,0,1],[0,0,1],[0,0,1],[0.01,0,1],[0.5,0.5,-1.0]],dtype=np.float64)
	
	hit_faces_control = np.array( [ 0,0, 1,1, -1 ], dtype=np.int64 )
	
	return vertices, faces, orgs, dirs, hit_faces_control



def test_TrimeshTracer():

	vertices, faces, orgs, dirs, hit_faces_control = get_test_data()
	
	tracer = POSBR.RayTracerTrimesh( vertices, faces )
	hit_faces = tracer.intersects_first( orgs, dirs )
	
	assert np.all( hit_faces == hit_faces_control )


def test_NanortTracer():

	vertices, faces, orgs, dirs, hit_faces_control = get_test_data()
	
	tracer = POSBR.RayTracerNanoRT( vertices, faces )
	hit_faces = tracer.intersects_first( orgs, dirs )
	
	assert np.all( hit_faces == hit_faces_control )

