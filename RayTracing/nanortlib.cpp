

#include "TracerNanoRT.hpp"
#include "Ray.hpp"
#include "Iterators.hpp"

#include<algorithm>
#include<execution>

extern "C" {

void* get_handle( unsigned long n_faces, double* vertices, unsigned* faces ){
	return (void*) new TracerNanoRT( n_faces, vertices, faces );
}

void intersects_first_interface( unsigned long n_rays, double* ray_orgs, double* ray_dirs, long* hit_indices, void* handle ){

	TracerNanoRT* tracer = (TracerNanoRT*)handle;
	
	// allocate space for t_fars
	std::vector<double> t_fars( n_rays, 1.0e+30 );
	// Wrap origins and directions
	Iterator2D org_it( ray_orgs, 3 );
	Iterator2D dir_it( ray_dirs, 3 );

	auto rays = RayListInterface( n_rays, org_it, dir_it, t_fars.begin(), hit_indices );
	
	std::for_each( std::execution::par, rays.begin(), rays.end(), *tracer );
}

void free_handle( void* handle ){
	TracerNanoRT* tracer = (TracerNanoRT*)handle;
	delete tracer;
}

} // extern "C"
