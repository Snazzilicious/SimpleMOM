
#ifndef RAY_HPP
#define RAY_HPP


/*
Specification for a "ray" to be traced by the Engine. Basic needs defined by tracer backend
Have not specified the "payload" of the ray (i.e. the info necessary to compute the ray-mesh interaction)
in case that is not universal.
Observations which depend on a certain type of payload should define their own ray which is castable
to & from this ray.
*/
template<typename OriginIterator, typename DirectionIterator>
class RayInterface {
	private:
		OriginIterator org_it;
		DirectionIterator dir_it;
		double* t_far_ptr;
		long* hit_face_ptr;
	
	public:
		RayInterface( const OriginIterator& org_begin, const DirectionIterator& dir_begin, double& tfar, long& hit_face ) : 
			org_it(org_begin), dir_it(dir_begin), t_far_ptr(&tfar), hit_face_ptr(&hit_face) {}
		
		RayInterface( const OriginIterator& org_begin, const DirectionIterator& dir_begin, double* tfar, long* hit_face ) : 
			org_it(org_begin), dir_it(dir_begin), t_far_ptr(tfar), hit_face_ptr(hit_face) {}
		
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
		RayListInterface( std::size_t n_rays, 
		                  OriginListIterator orgs_begin, 
		                  DirectionListIterator dirs_begin, 
		                  DoubleIterator t_fars_begin, 
		                  LongIterator hit_faces_begin )
		: _size( n_rays ), orgs( orgs_begin ), dirs( dirs_begin ), t_fars( t_fars_begin ), hit_faces( hit_faces_begin )
		{}
		
		std::size_t size(){ return _size; }
		
		struct iterator;
		
		iterator begin(){ return iterator( this ); }
		iterator end(){ return begin()+size(); }
		
		RayInterface<typename OriginListIterator::reference, typename DirectionListIterator::reference>
		operator[]( std::size_t i ){
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
		using value_type = RayInterface<typename OriginListIterator::reference, typename DirectionListIterator::reference> ;
		using pointer = RayInterface<typename OriginListIterator::reference, typename DirectionListIterator::reference> ;
		using reference = RayInterface<typename OriginListIterator::reference, typename DirectionListIterator::reference> ;
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


#endif /* RAY_HPP */

