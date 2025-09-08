

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








