
// TODO
// Design the leaves
// Constructing hMatrices
//    Assigning matrix data to Dense and Low rank blocks
//    Shared pointers to nodes + shared_from_this
//    How to create on the stack?
// How to know when pointing to leaf node? extra hMatrix node? add "type" labels? - leaning towards both
//    types for OOC and distributed matrices need to all be mutually disjoint, or need to be different node types
// Set item
//    need to identify and replace leaf nodes
//    LR should check if has same dataset, maybe
// subroutine error checking, debug and release, and early exiting
// Matrix multiply
//    All combinations
//        LR+LR should still check for same dataset
//    rSVD
// Tri Solves
//    should do the wrapper thing whenever encountering a LR matrix, then only solve when both dense
// Multithreaded
// OOC
//    slices must stop at the in-core block level i.e. a loadable hMatrix is the leaf data of an OOC hMatrix
//        loading a slice needs to be aware that its not pointing to the whole matrix
//    Ref counting on disk?
// Distributed Memory
//    how to "send" slices
// Distributed Disk
// API organization


void augment_low_rank( hMatrixInterface M, std::size_t new_rank ){
	// TODO
}

template<hMatrixSlice>
struct GEMM_job_descriptor{
	enum {GEMM,RSVD} job_type;
	hMatrixSlice a;
	hMatrixSlice b;
	hMatrixSlice c;
	GEMM_job_descriptor( t, a );
	GEMM_job_descriptor( t, a, b, c );
}

GEMM_job_descriptor queue_LR_H_LR_GEMM( hMatrixInterface A, hMatrixInterface B, hMatrixInterface C ){
	// augment C to receive result
	std::size_t old_rank = C.rank();
	std::size_t new_rank = old_rank;
	if( A.block_type() == LowRankMatrix )
		new_rank += A.rank();
	else
		new_rank += B.rank();
	augment_low_rank( C, new_rank );
	
	
	hMatrix c_wrapper, aorb_wrapper;
	hMatrixSlice a,b,c;
	if( A.block_type() == LowRankMatrix ){
		// TODO copy over idle basis
		C.left[:,old_rank:] = A.left[:,:] ;
		
		// TODO wrap remaining bases in hMatrices and pass along as operands
		aorb_wrapper.wrap( A.right[:,:] );
		c_wrapper.wrap( C.right[old_rank:,:] );
		
		a = aorb_wrapper.slice( 0,aorb_wrapper.nrows(), 0,aorb_wrapper.ncols() );
		b = B;
	}
	else {
		// copy over idle basis
		C.right[old_rank:,:] = B.right[:,:] ;
		
		// wrap remaining bases in hMatrices and pass along as operands
		aorb_wrapper.wrap( B.left[:,:] );
		c_wrapper.wrap( C.left[:,old_rank:] );
		
		a = A;
		b = aorb_wrapper.slice( 0,aorb_wrapper.nrows(), 0,aorb_wrapper.ncols() );
	}
	c = c_wrapper.slice( 0,c_wrapper.nrows(), 0,c_wrapper.ncols() );
	
	return GEMM_job_descriptor( GEMM, a, b, c );
}


void Leaf_GEMM( hMatrixInterface A, hMatrixInterface B, hMatrixInterface C ){
	// TODO
	auto type_a = a.block_type();
	auto type_b = b.block_type();
	auto type_c = c.block_type();
	
	if( type_a == DenseMatrix && type_b == DenseMatrix )
		// convert C to Dense
	
	if( type_c == DenseMatrix ){
		// going straight in
	
	}
	else {
		// augment C to receive result
	
		
	}
}


std::vector<std::size_t> get_partition_union( const std::vector<std::size_t>& part1, const std::vector<std::size_t>& part2 ){
	std::vector<std::size_t> part_union;
	std::set_union( part1.begin(), part1.end(), part2.begin(), part2.end(), std::back_inserter(part_union) );
	return part_union;
}

std::list<GEMM_job_descriptor> queue_H_GEMM( hMatrixInterface A, hMatrixInterface B, hMatrixInterface C ){
	// Get unions of partition points
	std::vector<std::size_t> row_begins = get_partition_union( A.get_row_begins(), C.get_row_begins() );
	std::vector<std::size_t> col_begins = get_partition_union( B.get_col_begins(), C.get_col_begins() );
	std::vector<std::size_t> inr_begins = get_partition_union( A.get_col_begins(), B.get_row_begins() );
	
	// Symbolically matmul all sub blocks
	std::list<GEMM_job_descriptor> job_stack;
	for( auto rb_it=row_begins.begin(); rb_it !=row_begins.end()-1; ++rb_it )
		for( auto cb_it=col_begins.begin(); cb_it !=col_begins.end()-1; ++cb_it )
			for( auto ib_it=inr_begins.begin(); ib_it !=inr_begins.end()-1; ++ib_it ){
	
				std::size_t row_begin = *rb_it;
				std::size_t row_end = *(rb_it+1);
				std::size_t col_begin = *cb_it;
				std::size_t col_end = *(cb_it+1);
				std::size_t inr_begin = *ib_it;
				std::size_t inr_end = *(ib_it+1);
				
				auto a = A.slice( row_begin, row_end, inr_begin, inr_end );
				auto b = B.slice( inr_begin, inr_end, col_begin, col_end );
				auto c = C.slice( row_begin, row_end, col_begin, col_end );
				// Partition C if needed to avoid slice invalidation by future operations
				if( c.root == C.root ){
					C.assign( row_begin, row_end, col_begin, col_end, c );
					c = C.slice( row_begin, row_end, col_begin, col_end );
				}
			
				job_stack.emplace_back( GEMM, a, b, c );
			}
		}
	}

	return job_stack;
}


template<Scalar,hMatrixInterface>
void hMatrixGEMM( Scalar alpha, hMatrixInterface A, hMatrixInterface B, hMatrixInterface C ){

	if ( !(A.nrows() == C.nrows() && B.ncols() == C.ncols() && A.ncols() == B.nrows()) )
		throw std::runtime_error("hMatrixGEMM: Incompatible operand dimensions.");

	if ( alpha == Scalar(0.0) ) return;
	
	std::list<GEMM_job_descriptor> job_stack;
	job_stack.emplace_back( GEMM, A, B, C );
	
	for( ; !job_stack.empty(); job_stack.pop_front() ){
	
		auto job = job_stack.front();
		
		if( job.type == RSVD ){
			rSVD( job.a );
		}
		else /* job.type == GEMM */ {
			
			auto a = job.a;
			auto b = job.b;
			auto c = job.c;
			
			assert( a.nrows() == c.nrows() && b.ncols() == c.ncols() && a.ncols() == b.nrows() );
			
			auto type_a = a.block_type();
			auto type_b = b.block_type();
			auto type_c = c.block_type();
			
			if( type_a != ZeroMatrix && type_b != ZeroMatrix ){
			
				// Special case of LowRankMatrix times hMatrix assigned to LowRankMatrix or ZeroMatrix
				if( type_c == ZeroMatrix || type_c == LowRankMatrix ){
					if( ( type_a == LowRankMatrix && type_b == hMatrix ) || ( type_a == hMatrix && type_b == LowRankMatrix ) ){
						GEMM_job_descriptor new_gemm = queue_LR_H_LR_GEMM( a, b, c );
						GEMM_job_descriptor new_rsvd = GEMM_job_descriptor( RSVD, c );
						// Push new jobs onto stack, right behind 'front'
						auto insert_pos = job_stack.begin();
						insert_pos++;
						insert_pos = job_stack.insert( insert_pos, new_rsvd );
						insert_pos = job_stack.insert( insert_pos, new_gemm );
					}
				}
				// Base case
				else if( type_a != hMatrix && type_b != hMatrix && type_c != hMatrix ){
					Leaf_GEMM( alpha, a, b, c );
				}
				// General case
				else {
					std::list<GEMM_job_descriptor> new_jobs = queue_H_GEMM( a, b, c );
					// Push new jobs onto stack, right behind 'front'
					auto insert_pos = job_stack.begin();
					insert_pos++;
					job_stack.insert( insert_pos, new_jobs.begin(), new_jobs.end() );
				}
			}
		}
	}
}



template<Scalar,hMatrixInterface>
void hMatrixGEMM( Scalar alpha, hMatrixInterface A, hMatrixInterface B, hMatrixInterface C, std::size_t n_threads ){

	// Split C up as fine as possible/reasonable
	std::list<GEMM_job_descriptor> job_stack;
	job_stack.emplace_back( GEMM, A, B, C );
	bool all_leaves = false;
	while( job_stack.size() < 2*n_threads && !all_leaves ){
		// Iterate through current queue and subdivide any jobs that can be
		all_leaves = true; 
		for( auto pos=job_stack.begin(); pos != job_stack.end(); ++pos ){
			
			auto a = pos->a;
			auto b = pos->b;
			auto c = pos->c;
			
			// if c can be subdivided into independent sections ...
			if( c.block_type() == hMatrix ){
				all_leaves = false;
				
				std::vector<std::size_t> row_begins = c.get_row_begins();
				std::vector<std::size_t> col_begins = c.get_col_begins();
				
				// Subdivide this job into new smaller ones
				std::list<GEMM_job_descriptor> new_jobs;
				for( auto rb_it=row_begins.begin(); rb_it !=row_begins.end()-1; ++rb_it ){
					for( auto cb_it=col_begins.begin(); cb_it !=col_begins.end()-1; ++cb_it ){
						
						std::size_t row_begin = *rb_it;
						std::size_t row_end = *(rb_it+1);
						std::size_t col_begin = *cb_it;
						std::size_t col_end = *(cb_it+1);
						
						auto new_a = a.slice( row_begin, row_end, 0, a.ncols() );
						auto new_b = b.slice( 0, b.nrows(), col_begin, col_end );
						auto new_c = C.slice( row_begin, row_end, col_begin, col_end );
						
						new_jobs.emplace_back( GEMM, new_a, new_b, new_c );
					}
				}
				
				// Replace current job with sub jobs
				pos = job_stack.erase( pos );
				job_stack.insert( pos, new_jobs.begin(), new_jobs.end() );
				--pos;
			}
		}
	}

	// Execute all operations in the queue concurrently
	#pragma omp parallel num_threads(n_threads) shared(job_stack)
	while( !job_stack.empty() ){
		bool do_job = false;
		GEMM_job_descriptor job;
		
		#pragma omp critical
		if( !job_stack.empty() ){
			job = job_stack.front();
			job_stack.pop_front();
			do_job = true;
		}
		
		if (do_job) hMatrixGEMM( alpha, job.a, job.b, job.c );
	}
}


template<Scalar,OOChMatrixInterface>
void hMatrixGEMM( Scalar alpha, OOChMatrixInterface A, OOChMatrixInterface B, OOChMatrixInterface C, std::size_t n_threads ){
	
	std::list<GEMM_job_descriptor> job_stack;
	job_stack.emplace_back( GEMM, A, B, C );
	
	// find loadable blocks
	for( ; !job_stack.empty(); job_stack.pop_front() ){
		
		auto job = job_stack.front();
		
		auto a = job.a;
		auto b = job.b;
		auto c = job.c;
		
		auto type_a = a.block_type();
		auto type_b = b.block_type();
		auto type_c = c.block_type();
		
		if( type_a != ZeroMatrix && type_b != ZeroMatrix ){
			// Base case
			if( type_a == IC_Block && type_b == IC_Block && type_c == IC_Block ){
				
				// load blocks - load so that a,b free data when they go out of scope
				auto incore_a = a.load();
				auto incore_b = b.load();
				auto incore_c = c.load();
				
				// do in-core GEMM
				hMatrixGEMM( alpha, incore_a, incore_b, incore_c, n_threads );
				
				// deallocate a,b / let go out of scope
				// TODO write c - c[:,:] = incore_c ???
			}
			// General case
			else {
				std::list<GEMM_job_descriptor> new_jobs = queue_H_GEMM( a, b, c );
				// Push new jobs onto stack, right behind 'front'
				auto insert_pos = job_stack.begin();
				insert_pos++;
				job_stack.insert( insert_pos, new_jobs.begin(), new_jobs.end() );
			}
		}
	}
}



template<Scalar,OOChMatrixInterface>
void hMatrixGEMM( Scalar alpha, OOChMatrixInterface A, OOChMatrixInterface B, OOChMatrixInterface C, std::size_t n_threads, MPI_comm Work_Comm ){

	// Chunk up C to find independent work
	
	// communicate tasking to each rank
	MPI_scatter( ... );
	
	// do OOC GEMM
	hMatrixGEMM( alpha, a, b, c, n_threads );
}


template<Scalar,OOChMatrixInterface>
void hMatrixGEMM( Scalar alpha, OOChMatrixInterface A, OOChMatrixInterface B, OOChMatrixInterface C,
                  std::size_t n_threads, MPI_comm Work_Comm, std::unordered_map<std::string,std::vector<int>> proc_groups ){

	// TODO
}

