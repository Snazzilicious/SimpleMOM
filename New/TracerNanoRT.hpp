
#ifndef TRACER_NANO_RT_HPP
#define TRACER_NANO_RT_HPP

#include "../External/nanort.h"

class TracerNanoRT {
	private:
		nanort::BVHAccel<double> accel;
		nanort::BVHTraceOptions trace_options; // default
		nanort::TriangleIntersector<double> triangle_intersector;
		
	public:
		TracerNanoRT( std::size_t n_facets, double* vertices_begin, unsigned* facets_begin ) : triangle_intersector( vertices_begin, facets_begin, sizeof(double) * 3 ) {
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


#endif /* TRACER_NANO_RT_HPP */
