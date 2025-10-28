
// TODO
// Tri Solves
//    should do the wrapper thing whenever encountering a LR matrix, then only solve when both dense


void Leaf_TRSM( char side, char uplo, char diag, hMatrixInterface A, hMatrixInterface B ){

	if ( !(A.block_type() == LU_Block && B.block_type() == DenseMatrix) )
		throw std::runtime_error("hMatrixTRSM: Invalid operand types.");
	
	// Get pointers
	
	// solve
	trsm( ... );
}

void hMatrixTRSM( char side, char uplo, char diag, hMatrixInterface A, hMatrixInterface B ){

	std::list<TRSM_job_descriptor> job_stack;
	job_stack.emplace_back( TRSM, A, B );
	
	for( ; !job_stack.empty(); job_stack.pop_front() ){
		
		auto job = job_stack.front();
		
		if( job.type == GEMM ){
			hMatrixGEMM( -1.0, job.a, job.b, job.c );
		}
		else /* job_type == TRSM */ {
		
			auto a = job.a;
			auto b = job.b;
			
			auto type_a = a.block_type();
			auto type_b = b.block_type();
			
			if( type_a == DenseMatrix && type_a == DenseMatrix ){
				Leaf_TRSM( ..., job.a, job.b );
			}
			else if ( type_b == LowRankMatrix ){
				// TODO wrap basis and continue
			}
			else {
				std::vector<std::size_t> diag_begins = get_partition_union( a.get_col_begins(), b.get_row_begins() );
				std::vector<std::size_t> col_begins = b.get_col_begins();
				
				auto insert_pos = job_stack.begin();
				insert_pos++;
				
				auto trsm_a = a.slice( diag_begins[0], diag_begins[1], diag_begins[0], diag_begins[1] );
				for( auto cb_it=col_begins.begin(); cb_it !=col_begins.end()-1; ++cb_it ){
					// create solve descriptors
					std::size_t col_begin = *cb_it;
					std::size_t col_end = *(cb_it+1);
					
					auto trsm_b = b.slice( diag_begins[0], diag_begins[1], col_begin, col_end );
					
					// push onto stack
					insert_pos = job_stack.insert( insert_pos, TRSM_job_descriptor( TRSM, trsm_a, trsm_b ) );
				}
				
				// create schur update descriptor
				auto gemm_a = a.slice( diag_begins[1], diag_begins.back(), 0, diag_begins[1] );
				auto gemm_b = b.slice( 0, diag_begins[1], 0, B.ncols() );
				auto gemm_c = c.slice( diag_begins[1], diag_begins.back(), 0, B.ncols() );
				insert_pos = job_stack.insert( insert_pos, TRSM_job_descriptor( GEMM, gemm_a, gemm_b, gemm_c ) );
				
				trsm_a = a.slice( diag_begins[1], diag_begins.back(), diag_begins[1], diag_begins.back() );
		
				insert_pos = job_stack.insert( insert_pos, TRSM_job_descriptor( TRSM, trsm_a, gemm_c ) );
			}
		}
	}
}


void hMatrixTRSM( char side, char uplo, char diag, hMatrixInterface A, hMatrixInterface B, std::size_t n_threads ){

	std::list<TRSM_job_descriptor> job_stack;
	job_stack.emplace_back( TRSM, A, B );
	
	for( ; !job_stack.empty(); job_stack.pop_front() ){
		
		auto job = job_stack.front();
		
		if( job.type == GEMM ){
			hMatrixGEMM( -1.0, job.a, job.b, job.c, n_threads );
		}
		else /* job_type == TRSM */ {
		
			auto a = job.a;
			auto b = job.b;
			
			auto type_a = a.block_type();
			auto type_b = b.block_type();
			
			if( type_a == DenseMatrix && type_a == DenseMatrix ){
				Leaf_TRSM( ..., a, b );
			}
			else if ( type_b == LowRankMatrix ){
				// TODO wrap basis and continue
			}
			else {
				std::vector<std::size_t> diag_begins = get_partition_union( a.get_col_begins(), b.get_row_begins() );
				std::vector<std::size_t> col_begins = b.get_col_begins();
				
				if( col_begins.size() > 0.75*n_threads ){
					#pragma omp parallel for num_threads(n_threads)
					for( auto cb_it=col_begins.begin(); cb_it !=col_begins.end()-1; ++cb_it ){
						std::size_t col_begin = *cb_it;
						std::size_t col_end = *(cb_it+1);
						
						auto rhs = b.slice( 0, b.nrows(), col_begin, col_end );
						
						hMatrixTRSM( ..., a, rhs );
					}
				}
				else{
					auto insert_pos = job_stack.begin();
					insert_pos++;
					
					auto trsm_a = a.slice( diag_begins[0], diag_begins[1], diag_begins[0], diag_begins[1] );
					for( auto cb_it=col_begins.begin(); cb_it !=col_begins.end()-1; ++cb_it ){
						// create solve descriptors
						std::size_t col_begin = *cb_it;
						std::size_t col_end = *(cb_it+1);
						
						auto trsm_b = b.slice( diag_begins[0], diag_begins[1], col_begin, col_end );
						
						// push onto stack
						insert_pos = job_stack.insert( insert_pos, TRSM_job_descriptor( TRSM, trsm_a, trsm_b ) );
					}
					
					// create schur update descriptor
					auto gemm_a = a.slice( diag_begins[1], diag_begins.back(), 0, diag_begins[1] );
					auto gemm_b = b.slice( 0, diag_begins[1], 0, B.ncols() );
					auto gemm_c = c.slice( diag_begins[1], diag_begins.back(), 0, B.ncols() );
					insert_pos = job_stack.insert( insert_pos, TRSM_job_descriptor( GEMM, gemm_a, gemm_b, gemm_c ) );
					
					trsm_a = a.slice( diag_begins[1], diag_begins.back(), diag_begins[1], diag_begins.back() );
			
					insert_pos = job_stack.insert( insert_pos, TRSM_job_descriptor( TRSM, trsm_a, gemm_c ) );
				}
			}
		}
	}
}
