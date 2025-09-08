
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



// Sample implementations of the above interfaces

#include<array>

struct Ray {
	protected:
		// Basic needs
		std::array<double,3> origin;
		std::array<double,3> direction;
		double _t_far;
		long hit_face;
	
	public:
		Ray() : hit_face(-1), _t_far(1.0e30) {}
		Ray( std::array<double,3> org, std::array<double,3> dir ) : Ray() {
			origin = org;
			direction = dir;
		}
		
		// Core interface
		double* origin_begin(){ return origin.data(); }
		double* direction_begin(){ return direction.data(); }
		double& t_far(){ return _t_far; }
		long& hit_face_ind(){ return hit_face; }
		
		const double* origin_begin() const { return origin.data(); }
		const double* direction_begin() const { return direction.data(); }
		const double& t_far() const { return _t_far; }
		const long& hit_face_ind() const { return hit_face; }
		
		template<typename ray_t>
		void operator=( ray_t r_in ){
			std::copy( r_in.origin_begin(),r_in.origin_begin()+3, origin_begin() );
			std::copy( r_in.direction_begin(),r_in.direction_begin()+3, direction_begin() );
			t_far() = r_in.t_far();
			hit_face_ind() = r_in.hit_face_ind();
		}
};

#include<vector>

class RayList {
	private:
		std::vector<double> origins;
		std::vector<double> directions;
		std::vector<double> t_fars;
		std::vector<long> hit_faces;

	public:
		RayList( std::size_t length ){ resize(length); }	

		std::size_t size(){ return t_fars.size(); }
		
		void resize( std::size_t new_length ){
			origins.resize( 3*new_length );
			directions.resize( 3*new_length );
			t_fars.resize( new_length );
			hit_faces.resize( new_length );
		}

		struct iterator;
		
		iterator begin();
		iterator end();
		
		friend struct RayRef;
};


#include "Iterators.hpp"

struct RayRef {
	private:
		RayList* container;
		std::size_t ind;
	
	public:
		RayRef( RayList* parent, std::size_t index=0 ) : container(parent), ind(index) {}
		
		StrideIterator<double*> origin_begin(){ return StrideIterator<double*>( container->origins.data()+ind, container->size() ); }
		StrideIterator<double*> direction_begin(){ return StrideIterator<double*>( container->directions.data()+ind, container->size() ); }
		double& t_far(){ return container->t_fars[ind]; }
		long& hit_face_ind(){ return container->hit_faces[ind]; }
		
		template<typename ray_t>
		void operator=( ray_t r_in ){
			std::copy( r_in.origin_begin(),r_in.origin_begin()+3, origin_begin() );
			std::copy( r_in.direction_begin(),r_in.direction_begin()+3, direction_begin() );
			t_far() = r_in.t_far();
			hit_face_ind() = r_in.hit_face_ind();
		}
};


struct RayList::iterator {
	private:
		RayList* container;
		std::size_t pos;
		
	public:
		using difference_type = std::ptrdiff_t ;
		using value_type = Ray ;
		using pointer = RayRef ;
		using reference = RayRef ;
		using iterator_category = std::random_access_iterator_tag ;

		iterator( RayList* parent, std::size_t index=0 ) : container(parent), pos(index) {}
		
		iterator operator++(){ ++pos; return *this; }
		iterator operator+(difference_type n){ iterator tmp = *this; tmp.pos += n; return tmp; }
		iterator& operator+=(difference_type n){ pos += n; return *this; }
		friend difference_type operator-(const iterator& a, const iterator& b){ return a.pos - b.pos; }
		
		friend bool operator==( const iterator& a, const iterator& b ){ return a.pos == b.pos; }
		friend bool operator!=( const iterator& a, const iterator& b ){ return !(a == b); }
		
		reference operator[](difference_type n){ return *(operator+(n)); }
		reference operator*(){ return RayRef( container, pos ); }
		pointer operator->(){ return RayRef( container, pos ); }
};


#endif /* RAY_HPP */

