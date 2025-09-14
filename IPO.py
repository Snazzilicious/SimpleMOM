
import numpy as np


from . import POSBR

def is_illuminated_pw( prop_dir, normals, tol=1e-1 ):
	return np.sum( normals * prop_dir, axis=1 ) <= -tol


def IPO_current_from_H( prop_dir, H_pol, normals, tol=1e-1 ):
	illuminated = is_illuminated_pw( prop_dir, normals, tol )
	currents = POSBR.PO_current_from_H( H_pol, normals )
	return ( currents.T * illuminated ).T


def IPO_current_from_E( prop_dir, E_pol, normals, tol=1e-1 ):
	illuminated = is_illuminated_pw( prop_dir, normals, tol )
	currents = POSBR.PO_current_from_E( prop_dir, E_pol, normals )
	return ( currents.T * illuminated ).T




""" Make cluster tree """
class InnerNode:
	def __init__( self, begin_ind, end_ind, sorted_inds, pts ):
		self.children = None
		
		self.begin = begin_ind
		self.end = end_ind
		
		self.center = np.mean( pts[sorted_inds[self.begin:self.end]], axis=0 )
		self.R = np.max( np.linalg.norm( pts[sorted_inds[self.begin:self.end]] - self.center, axis=1 ) )
		
	def n_members(self):
		return self.end - self.begin


from sklearn.cluster import KMeans

def get_cluster_tree( centroids ):
	N_CHILDREN = 20
	MAX_MEMBERS = 100

	sorted_inds = np.arange(len(centroids))

	node_stack = [ InnerNode( 0, len(centroids), sorted_inds, centroids ) ]
	root = node_stack[0]

	while len(node_stack) > 0:
		node = node_stack.pop()
		
		# Group spatially
		km = KMeans( n_clusters=N_CHILDREN, init='k-means++', n_init=1 ).fit(centroids[sorted_inds[node.begin:node.end]])
		pivot = np.argsort( km.labels_ )
		ind_ranges = [node.begin]
		ind_ranges.extend([ node.begin + i for i in np.cumsum( np.bincount(km.labels_) ) ])
		# 
		sorted_inds[node.begin:node.end] = sorted_inds[node.begin:node.end][pivot]

		# make new nodes for each group
		node.children = [ InnerNode( ind_ranges[i], ind_ranges[i+1], sorted_inds, centroids ) for i in range(N_CHILDREN) ]

		# push all the children onto a stack/queue
		node_stack.extend([ child for child in node.children if child.n_members() >= MAX_MEMBERS ])

	return root, sorted_inds



""" Construct IPO visibility matrix (using cluster tree) """
def get_recieve_map( centroids, normals, tol=1e-1 ):
	root, sorted_inds = get_cluster_tree( centroids )
	
	recv_map = [[] for _ in range(len(centroids))]

	for face_ind in range(len(centroids)):
		node_stack = [root]
		
		while len(node_stack) > 0:
			node = node_stack.pop()
			
			if node.children is None: # terminal node
				centers = centroids[sorted_inds[node.begin:node.end]]
				vs = centroids[face_ind] - centers
				norm_vs = np.linalg.norm( vs, axis=1 )
				nTvs = vs.dot(normals[face_ind]) / (norm_vs + 1e-12)
				# face recieves from those with sufficiently negative (normalized) dot product
				recv_map[face_ind].extend( sorted_inds[node.begin:node.end][nTvs <= -tol] )
				
			else:
				for child in node.children:
					v = centroids[face_ind] - child.center
					norm_v = np.linalg.norm(v)
					
					if child.R >= norm_v :
						# must split further
						node_stack.append(child)
					else:
						nTv = normals[face_ind].dot( v )
						tTv = np.sqrt( norm_v**2 - nTv**2 )
						a_p = ( child.R / norm_v**2 ) * ( -child.R*nTv + tTv*np.sqrt( norm_v**2 - child.R**2 ) )
						a_m = ( child.R / norm_v**2 ) * ( -child.R*nTv - tTv*np.sqrt( norm_v**2 - child.R**2 ) )
						b_p = -( a_p*nTv + child.R**2 ) / tTv
						b_m = -( a_m*nTv + child.R**2 ) / tTv
						L_p = ( nTv + a_p ) / np.sqrt( norm_v**2 + child.R**2 + 2*( a_p*nTv + b_p*tTv ) )
						L_m = ( nTv + a_m ) / np.sqrt( norm_v**2 + child.R**2 + 2*( a_m*nTv + b_m*tTv ) )
						
						L_max = max(L_p,L_m)
						L_min = min(L_p,L_m)
						if L_max <= -tol and L_min <= -tol :
							# face recieves from all in group
							recv_map[face_ind].extend( sorted_inds[child.begin:child.end] )
						elif L_max > -tol and L_min > -tol :
							# face recieves from none
							pass
						else:
							# indeterminate - must split further
							node_stack.append(child)
		
		recv_map[face_ind] = np.array( recv_map[face_ind], dtype=np.int64 )
		recv_map[face_ind].sort()

	return recv_map


def get_recieve_map_optimized( centroids, normals, tol=1e-1 ):
	root, sorted_inds = get_cluster_tree( centroids )
	
	recv_map = [[] for _ in range(len(centroids))]

	node_stack = [(root,np.arange(len(centroids)))]

	while len(node_stack) > 0:
		node,inds = node_stack.pop()
		
		if node.children is None: # terminal node
			centers = centroids[sorted_inds[node.begin:node.end]]
			for face_ind in inds:
				vs = centroids[face_ind] - centers
				norm_vs = np.linalg.norm( vs, axis=1 )
				nTvs = vs.dot(normals[face_ind]) / (norm_vs+1e-12)
				recv_map[face_ind].extend( sorted_inds[node.begin:node.end][nTvs <= -tol] )
		
		else:
			for child in node.children:
				v = centroids[inds] - child.center
				norm_v = np.linalg.norm( v, axis=1 )
				
				inside = child.R >= norm_v
				outside = np.logical_not(inside)
				descend_inds = inds[inside].tolist() # must split further
				check_inds = inds[outside]
				
				v = v[outside]
				norm_v = norm_v[outside]
				
				nTv = np.sum( normals[check_inds] * v, axis=1 )
				tTv = np.sqrt( norm_v**2 - nTv**2 )
				a_p = ( child.R / norm_v**2 ) * ( -child.R*nTv + tTv*np.sqrt( norm_v**2 - child.R**2 ) )
				a_m = ( child.R / norm_v**2 ) * ( -child.R*nTv - tTv*np.sqrt( norm_v**2 - child.R**2 ) )
				b_p = -( a_p*nTv + child.R**2 ) / tTv
				b_m = -( a_m*nTv + child.R**2 ) / tTv
				L_p = ( nTv + a_p ) / np.sqrt( norm_v**2 + child.R**2 + 2*( a_p*nTv + b_p*tTv ) )
				L_m = ( nTv + a_m ) / np.sqrt( norm_v**2 + child.R**2 + 2*( a_m*nTv + b_m*tTv ) )
				L_max = np.maximum( L_p, L_m )
				L_min = np.minimum( L_p, L_m )
				
				# faces recieves from all in group
				for face_ind in check_inds[ np.logical_and( L_max <= -tol, L_min <= -tol ) ]:
					recv_map[face_ind].extend( sorted_inds[child.begin:child.end] )
				
				# indeterminate - must split further
				descend_inds.extend( check_inds[ np.logical_xor( L_max <= -tol, L_min <= -tol ) ] )
				if len(descend_inds) > 0:
					node_stack.append( ( child, np.array(descend_inds) ) )
		
	for i in range(len(recv_map)):
		recv_map[i] = np.array( recv_map[i], dtype=np.int64 )
		recv_map[i].sort()
	
	return recv_map


def get_recieve_map_brute_force( centroids, normals, tol=1e-1 ):
	recv_map = []
	for face_ind in range(len(centroids)):
		vs = centroids[face_ind] - centroids
		norm_vs = np.linalg.norm( vs, axis=1 )
		nTvs = vs.dot(normals[face_ind]) / (norm_vs+1e-12)
		recv_map.append( np.where( nTvs <= -tol )[0] )
	
	return recv_map



""" Diagnostics """
def n_mismatch( map1, map2 ):
	s=0
	for a,b in zip(map1,map2):
		for i,j in zip(a,b):
		    s += i != j
	return s

def n_nonunique( map1 ):
	s=0
	for l in map1:
		s += np.any( l != np.unique(l) )
	return s

def get_clusters( root, sorted_inds ):
	node_stack=[root]
	cluster_labels=[]
	while len(node_stack) > 0:
		cluster_labels.append(-np.ones(len(sorted_inds)))
		new_node_stack = []
		for i,node in enumerate(node_stack):
			cluster_labels[-1][sorted_inds[node.begin:node.end]] = i
			if node.children is not None:
				new_node_stack.extend(node.children)
		
		node_stack=new_node_stack
	
	return cluster_labels


