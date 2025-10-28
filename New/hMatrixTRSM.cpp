
// TODO
// Tri Solves
//    should do the wrapper thing whenever encountering a LR matrix, then only solve when both dense
//    if low rank blocks share a basis, must not solve twice


void Leaf_TRSM( char side, char uplo, char diag, hMatrixInterface A, hMatrixInterface B ){

	if ( !(A.block_type() == LU_Block && B.block_type() == DenseMatrix) )
		throw std::runtime_error("hMatrixTRSM: Invalid operand types.");
	
	// Get pointers
	
	// solve
	trsm( ... );
}



std::list<TRSM_job_descriptor> queue_H_TRSM( hMatrixInterface A, hMatrixInterface B ){

	std::vector<std::size_t> diag_begins = get_partition_union( a.get_col_begins(), b.get_row_begins() );
	std::vector<std::size_t> col_begins = b.get_col_begins()
	
	for( auto db_it=diag_begins.begin(); db_it !=diag_begins.end()-1; ++db_it ){
		for( auto cb_it=col_begins.begin(); cb_it !=col_begins.end()-1; ++cb_it ){
			// create solve descriptors
			std::size_t diag_begin = *db_it;
			std::size_t diag_end = *(db_it+1);
			std::size_t col_begin = *cb_it;
			std::size_t col_end = *(cb_it+1);
			
			
			
		}
		
		// create schur update descriptor
		auto gemm_a = a.slice( diag_begins[1], diag_begins.back(), 0, diag_begins[1] );
		auto gemm_b = b.slice( 0, diag_begins[1], 0, B.ncols() );
		auto gemm_c = c.slice( diag_begins[1], diag_begins.back(), 0, B.ncols() );
		TRSM_job_descriptor new_gemm( GEMM, gemm_a, gemm_b, gemm_c );
	}

	// Get unions of partition points
	std::vector<std::size_t> row_begins = get_partition_union( A.get_row_begins(), C.get_row_begins() );
	std::vector<std::size_t> col_begins = get_partition_union( B.get_col_begins(), C.get_col_begins() );
	std::vector<std::size_t> inr_begins = get_partition_union( A.get_col_begins(), B.get_row_begins() );
	
	// Symbolically matmul all sub blocks
	std::list<GEMM_job_descriptor> job_stack;
	for( auto rb_it=row_begins.begin(); rb_it !=row_begins.end()-1; ++rb_it ){
		for( auto cb_it=col_begins.begin(); cb_it !=col_begins.end()-1; ++cb_it ){
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


void hMatrixTRSM( char side, char uplo, char diag, hMatrixInterface A, hMatrixInterface B ){

	std::list<TRSM_job_descriptor> job_stack;
	job_stack.emplace_back( TRSM, A, B );
	
	for( ; !job_stack.empty(); job_stack.pop_front() ){
		
		auto job = job_stack.front();
		
		if( job.type == GEMM ){
			hMatrixGEMM( -1.0, job.a, job.b, job.c );
		}
		else /* job_type == TRSM */ {
			
			if( both dense ){
				Leaf_TRSM( ..., job.a, job.b );
			}
			else if ( type_b == LowRankMatrix ){
				// wrap basis and continue
			}
			else {
				
				
				// push onto stack
				auto insert_pos = job_stack.begin();
				insert_pos++;
				insert_pos = job_stack.insert( insert_pos, new_gemm );
				insert_pos = job_stack.insert( insert_pos, new_trsm );
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
			
			if( both dense ){
				Leaf_TRSM( ..., job.a, job.b );
			}
			else if ( type_b == LowRankMatrix ){
				// wrap basis and continue
			}
			else {
				
				
				// push onto stack
				auto insert_pos = job_stack.begin();
				insert_pos++;
				insert_pos = job_stack.insert( insert_pos, new_gemm );
				insert_pos = job_stack.insert( insert_pos, new_trsm );
			}
		}
	}
}
