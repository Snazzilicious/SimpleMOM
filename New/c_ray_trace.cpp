


template<typename OriginIterator, typename DirectionIterator>
class RayInterface {
	private:
		OriginIterator org_it;
		DirectionIterator dir_it;
		double* t_far_ptr;
		long* hit_face_ptr;
	
	public:
		RayInterface( const OriginIterator& org_begin, const DirectionIterator& dir_begin, double& t_far, long& hit_face ) : 
			org_it(org_begin), dir_it(dir_begin), t_far_ptr(&t_far), hit_face_ptr(&hit_face) {}
		
		RayInterface( const OriginIterator& org_begin, const DirectionIterator& dir_begin, double* t_far, long* hit_face ) : 
			org_it(org_begin), dir_it(dir_begin), t_far_ptr(t_far), hit_face_ptr(hit_face) {}
		
		// Core ray interface
		OriginIterator origin_begin(){ return org_it; }
		DirectionIterator direction_begin(){ return dir_it; }
		double& t_far(){ return *t_far_ptr; }
		long& hit_face_ind(){ return *hit_face_ptr; }
		
		template<typename ray_t>
		void operator=( ray_t r_in ){
			std::copy( r_in.origin_begin(),r_in.origin_begin()+3, origin_begin() );
			std::copy( r_in.direction_begin(),r_in.direction_begin()+3, direction_begin() );
			t_far() = r_in.t_far();
			hit_face_ind() = r_in.hit_face_ind();
		}
};

template<typename OriginListIterator, typename DirectionListIterator, typename DoubleIterator, typename LongIterator>
class RayListInterface {
	private:
		OriginListIterator orgs;
		DirectionListIterator dirs;
		DoubleIterator t_fars;
		LongIterator hit_faces;
		std::size_t _size;

	public:
		
		std::size_t size(){ return _size; }
		
		struct iterator;
		
		iterator begin(){ return iterator( this ); }
		iterator end(){ return begin()+size(); }
		
		RayInterface operator[]( std::size_t i ){
			return RayInterface( orgs[i], dirs[i], t_fars[i], hit_faces[i] );
		}
};

template<typename OriginListIterator, typename DirectionListIterator, typename DoubleIterator, typename LongIterator>
struct RayListInterface<OriginListIterator, DirectionListIterator, DoubleIterator, LongIterator>::iterator {
	private:
		RayListInterface<OriginListIterator, DirectionListIterator, DoubleIterator, LongIterator>* container;
		std::size_t pos;
		
	public:
		using difference_type = std::ptrdiff_t ;
		using value_type = RayInterface ;
		using pointer = RayInterface ;
		using reference = RayInterface ;
		using iterator_category = std::random_access_iterator_tag ;

		iterator( RayListInterface<OriginListIterator, DirectionListIterator, DoubleIterator, LongIterator>* parent, std::size_t index=0 ) 
			: container(parent), pos(index) {}
		
		iterator operator++(){ ++pos; return *this; }
		iterator operator+(difference_type n){ iterator tmp = *this; tmp.pos += n; return tmp; }
		iterator& operator+=(difference_type n){ pos += n; return *this; }
		friend difference_type operator-(const iterator& a, const iterator& b){ return a.pos - b.pos; }
		
		friend bool operator==( const iterator& a, const iterator& b ){ return a.pos == b.pos; }
		friend bool operator!=( const iterator& a, const iterator& b ){ return !(a == b); }
		
		reference operator[](difference_type n){ return *(operator+(n)); }
		reference operator*(){ return container->operator[]( pos ); }
		pointer operator->(){ return operator*(); }
};




#include "../External/nanort.h"

class TracerNanoRT {
	private:
		nanort::BVHAccel<double> accel;
		nanort::BVHTraceOptions trace_options; // default
		nanort::TriangleIntersector<double> triangle_intersector;
		
	public:
		TracerNanoRT( std::size_t n_facets, double* vertices_begin, unsigned* facets_begin )
			: triangle_intersector( vertices_begin, facets_begin, sizeof(double) * 3 )
		{
			// Wrap the mesh in the required objects
			std::size_t vertices_stride = sizeof(double) * 3;
			nanort::TriangleMesh<double> triangle_mesh( vertices_begin, facets_begin, vertices_stride );
			nanort::TriangleSAHPred<double> triangle_pred( vertices_begin, facets_begin, vertices_stride );
			nanort::BVHBuildOptions<double> options; // default
			bool success = accel.Build( n_facets, triangle_mesh, triangle_pred, options );
			assert(success);
		}
		
		template<typename ray_t>
		void operator()(ray_t r_in){
			// Set ray in appropriate format
			nanort::Ray<double> ray;
			ray.min_t = 0.0;
			ray.max_t = 1.0e+30;
			ray.org[0] = r_in.origin_begin()[0];
			ray.org[1] = r_in.origin_begin()[1];
			ray.org[2] = r_in.origin_begin()[2];
			ray.dir[0] = r_in.direction_begin()[0];
			ray.dir[1] = r_in.direction_begin()[1];
			ray.dir[2] = r_in.direction_begin()[2];

			nanort::TriangleIntersection<double> isect;
			bool hit = accel.Traverse( ray, triangle_intersector, &isect, trace_options );
			if (hit) {
				r_in.hit_face_ind() = isect.prim_id;
				r_in.t_far() = isect.t;
			}
		}
};


// TODO Tracer type for python


#include "Iterators.hpp"

// struck first: index of struck or -1
void first_hit( std::size_t n_rays, double* origins, double* dirs ){
	
	// allocate tFars struck faces containers
	std::vector<double> t_fars( n_rays );
	std::vector<long> hit_faces( n_rays, -1 );

	// create appropiate wrapper for origins and dirs
	
	// build ray list
	RayListInterface ray_buffer( n_rays, origins_wrapper, dirs_wrapper, t_fars.begin(), hit_faces.begin() );
	
	// Trace rays
	std::for_each( ray_buffer.begin(), ray_buffer.end(), nanort_trace );
}

// TODO Python interface








