
#ifndef ENGINE_HPP
#define ENGINE_HPP

/*
Main class for tying together the basic building blocks mesh, excitations, observations etc.
*/


#include<mpi.h>

template<typename List>
std::size_t reduce_ray_queue_sizes( const List& queue, MPI_Comm comm ){
	unsigned long length = queue.size();
	unsigned long reduced_length;
	(void)MPI_Allreduce( &length, &reduced_length, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm );
	return reduced_length;
}


#include "../Common/Iterators.hpp"
#include "Ray.hpp"
#include "MeshInterface.hpp"
#include "TracerNanoRT.hpp"


namespace Engine {

template<class RayGenerator, class Observation>
void trace_excitation( const MeshInterface& mesh, RayGenerator& src, Observation& obs, int piece_ind=0, int n_pieces=1, MPI_Comm comm=MPI_COMM_WORLD )
{
	int proc_id, n_proc;
	(void)MPI_Comm_size( comm, &n_proc );
	(void)MPI_Comm_rank( comm, &proc_id );
	
	// Get this process's ray range
	// Allocate work associated with this piece
	auto piece_range = distribute_round_robin( src.rays_begin(), src.rays_end(), n_pieces, piece_ind );
	// and distrbute across processes
	auto proc_range = distribute_round_robin( piece_range.first, piece_range.second, n_proc, proc_id );
	StrideIterator<StrideIterator<typename RayGenerator::iterator>> src_begin = proc_range.first;
	StrideIterator<StrideIterator<typename RayGenerator::iterator>> src_end = proc_range.second;
	std::size_t n_src_rays = src_end - src_begin ;
	const std::size_t BATCH_SIZE = 1024 ;
	
	// Get initial batch of rays
	typename RayGenerator::RayList src_ray_buffer( BATCH_SIZE );
	std::copy_n( src_begin, BATCH_SIZE, src_ray_buffer.begin() );
	src_begin += BATCH_SIZE ;
	n_src_rays -= BATCH_SIZE ;
	
	
	// Set up Ray tracer (nanort)
	TracerNanoRT nanort_trace( mesh.n_facets(), mesh.vertices_begin()[0], mesh.facets_begin()[0] );
	
	while( reduce_ray_queue_sizes( src_ray_buffer, comm ) > 0 ){ // this check is MPI_Communicator wide
	
		// Format how tracer wants XXX if no significant reformatting, tracer loop could draw directly from src_ray_buffer XXX
		RayList tracer_ray_buffer( BATCH_SIZE );
		std::copy( src_ray_buffer.begin(), src_ray_buffer.end(), tracer_ray_buffer.begin() );
		
		// Trace rays
		std::for_each( tracer_ray_buffer.begin(), tracer_ray_buffer.end(), nanort_trace );
		
		// Write back to source format XXX if no significant reformatting, tracer loop could write directly to src_ray_buffer XXX
		std::copy( tracer_ray_buffer.begin(), tracer_ray_buffer.end(), src_ray_buffer.begin() );
		
		// Invoke hit handler
		obs.process_traced_rays( src_ray_buffer, mesh, comm );
		
		
		// Keep source ray buffer full
		if( src_ray_buffer.size() < BATCH_SIZE ){
			std::size_t n_rays = src_ray_buffer.size();
			src_ray_buffer.resize( std::min( BATCH_SIZE, n_rays+n_src_rays ) );
			std::size_t n_more_rays = src_ray_buffer.size()-n_rays;
			std::copy_n( src_begin, n_more_rays, src_ray_buffer.begin()+n_rays );
			src_begin += n_more_rays ;
			n_src_rays -= n_more_rays ;
		}
	}
	
	obs.global_reduce_all( comm );
}

} // End namespace Engine


#endif /* ENGINE_HPP */
