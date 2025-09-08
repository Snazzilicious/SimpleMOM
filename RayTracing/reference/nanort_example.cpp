

//#define NANORT_USE_CPP11_FEATURE // Enable C++11 feature
//#define NANORT_ENABLE_PARALLEL_BUILD // Enable parallel BVH build(OpenMP version is not yet fully tested)
#include "External/nanort.h"

#include<iostream>
#include<omp.h>

int main( int argc, char **argv ){

	// load mesh data...
	std::vector<float> vertices = {0.0,0.0,0.0, 1.0,0.2,0.0, 0.0,1.0,0.0};
	std::vector<unsigned int> facets = {0,1,2};
	
	// Wrap the mesh in required objects
	std::size_t stride = sizeof(float) * 3;
	nanort::TriangleMesh<float> triangle_mesh( vertices.data(), facets.data(), stride );
	nanort::TriangleSAHPred<float> triangle_pred( vertices.data(), facets.data(), stride );
	nanort::BVHBuildOptions<float> options; // Use default options
	nanort::BVHAccel<float> accel;
	bool success = accel.Build( facets.size(), triangle_mesh, triangle_pred, options );
	assert(success);

	
	const float tFar = 1.0e+30f;
	nanort::BVHTraceOptions trace_options;
	nanort::TriangleIntersector<> triangle_intersector( vertices.data(), facets.data(), stride );
	
	// Shoot rays.
//	#pragma omp parallel for
	for (int y = 0; y < 3; y++) {
		// Set ray in appropriate format
		nanort::Ray<float> ray;
		ray.min_t = 0.0f;
		ray.max_t = tFar;
		ray.org[0] = 0.4f;
		ray.org[1] = 0.4f;
		ray.org[2] = 1.0f;
		ray.dir[0] = -0.1;
		ray.dir[1] = -0.1;
		ray.dir[2] = -1.0 + 0.75*y;

		nanort::TriangleIntersection<> isect;
		bool hit = accel.Traverse( ray, triangle_intersector, &isect, trace_options );
		if (hit) {
			// Handle hit
			unsigned int fid = isect.prim_id;
			std::cout << omp_get_thread_num() << ": Struck face " << fid << std::endl;
			std::cout << "Local Coordinates (" << isect.u << ", " << isect.v << ")" << std::endl;
			std::cout << "Distance along dir: " << isect.t << std::endl;
		}
	}
	
	return 0;
}
