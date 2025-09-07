
#ifndef RAY_HPP
#define RAY_HPP

#include<array>

/*
Specification for a "ray" to be traced by the Engine. Basic needs defined by tracer backend
Have not specified the "payload" of the ray (i.e. the info necessary to compute the ray-mesh interaction)
in case that is not universal.
Observations which depend on a certain type of payload should define their own ray which is castable
to & from this ray.
*/
struct Ray {
	protected:
		// Basic needs
		std::array<double,3> origin;
		std::array<double,3> direction;
		double tFar;
		long hit_face;
	
	public:
		Ray() : hit_face(-1), tFar(1.0e30) {}
		Ray( std::array<double,3> org, std::array<double,3> dir ) : Ray() {
			origin = org;
			direction = dir;
		}
		
		// Core interface
		double* origin_begin(){ return origin.data(); }
		double* direction_begin(){ return direction.data(); }
		double& t_far(){ return tFar; }
		long& hit_face_ind(){ return hit_face; }
		
		const double* origin_begin() const { return origin.data(); }
		const double* direction_begin() const { return direction.data(); }
		const double& t_far() const { return tFar; }
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
		std::vector<double> tFars;
		std::vector<long> hit_faces;

	public:
		RayList( std::size_t length ){ resize(length); }	

		std::size_t size(){ return tFars.size(); }
		
		void resize( std::size_t new_length ){
			origins.resize( 3*new_length );
			directions.resize( 3*new_length );
			tFars.resize( new_length );
			hit_faces.resize( new_length );
		}

		struct iterator;
		
		iterator begin();
		iterator end();
		
		friend struct RayRef;
};


#include "../Common/Iterators.hpp"

struct RayRef {
	private:
		RayList* container;
		std::size_t ind;
	
	public:
		RayRef( RayList* parent, std::size_t index=0 ) : container(parent), ind(index) {}
		
		StrideIterator<double*> origin_begin(){ return StrideIterator<double*>( container->origins.data()+ind, container->size() ); }
		StrideIterator<double*> direction_begin(){ return StrideIterator<double*>( container->directions.data()+ind, container->size() ); }
		double& t_far(){ return container->tFars[ind]; }
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

