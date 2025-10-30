
// TODO
// Tri Solves
//    need ?laswp at some level
//        Probably needs to be at the hMatrix level


void Leaf_TRSM( char side, char uplo, char diag, hMatrixInterface A, hMatrixInterface B ){

	if ( !(A.block_type() == LU_Block && B.block_type() == DenseMatrix) )
		throw std::runtime_error("hMatrixTRSM: Invalid operand types.");
	
	// Get pointers of dense data, (and strides and layouts)
	
	// solve
	trsm( ... );
}


std::list<TRSM_job_descriptor> queue_H_TRSM( char side, char uplo, A, B ){

	std::list<TRSM_job_descriptor> job_stack;
	if( side == "L" && uplo == "L" ){
	
		std::vector<std::size_t> diag_begins = B.get_row_begins();
		std::vector<std::size_t> rhs_begins = B.get_col_begins();
			
		for( auto db_it=diag_begins.begin(); db_it != diag_begins.end()-1; ++db_it ){
			
			std::size_t diag_begin = *db_it;
			std::size_t diag_end = *(db_it+1);
			
			for( auto rb_it=rhs_begins.begin(); cb_it != rhs_begins.end()-1; ++cb_it ){
				
				std::size_t rhs_begin = *rb_it;
				std::size_t rhs_end = *(rb_it+1);
				
				auto a = A.slice( diag_begin, diag_end, diag_begin, diag_end );
				auto b = B.slice( diag_begin, diag_end, rhs_begin, rhs_end );
				
				job_stack.emplace_back( TRSM, a, b );
			}
			
			auto a = A.slice( diag_end, diag_begins.back(), diag_begin, diag_end );
			auto b = B.slice( diag_begin, diag_end, 0, rhs_begins.back() );
			auto c = B.slice( diag_end, diag_begins.back(), 0, rhs_begins.back() );
			
			job_stack.emplace_back( GEMM, a, b, c );
		}
	}
	else if( side == "R" && uplo == "U" ) {
	
		std::vector<std::size_t> diag_begins = B.get_col_begins();
		std::vector<std::size_t> rhs_begins = B.get_row_begins();
			
		for( auto db_it=diag_begins.begin(); db_it != diag_begins.end()-1; ++db_it ){
			
			std::size_t diag_begin = *db_it;
			std::size_t diag_end = *(db_it+1);
			
			for( auto rb_it=rhs_begins.begin(); cb_it != rhs_begins.end()-1; ++cb_it ){
				
				std::size_t rhs_begin = *rb_it;
				std::size_t rhs_end = *(rb_it+1);
				
				auto a = A.slice( diag_begin, diag_end, diag_begin, diag_end );
				auto b = B.slice( rhs_begin, rhs_end, diag_begin, diag_end );
				
				job_stack.emplace_back( TRSM, a, b );
			}
			
			auto a = B.slice( 0, rhs_begins.back(), diag_begin, diag_end );
			auto b = A.slice( diag_begin, diag_end, diag_end, diag_begins.back() );
			auto c = B.slice( 0, rhs_begins.back(), diag_end, diag_begins.back() );
			
			job_stack.emplace_back( GEMM, a, b, c );
		}
	}
	else if( side == "L" && uplo == "U" ){
		
		std::vector<std::size_t> diag_begins = B.get_row_begins();
		std::vector<std::size_t> rhs_begins = B.get_col_begins();
		
		for( auto db_it=diag_begins.end()-1; db_it != diag_begins.begin(); --db_it ){
			
			std::size_t diag_begin = *(db_it-1);
			std::size_t diag_end = *db_it;
			
			auto a = A.slice( diag_begin, diag_end, diag_end, diag_begins.back() );
			auto b = B.slice( diag_end, diag_begins.back(), 0, rhs_begins.back() );
			auto c = B.slice( diag_begin, diag_end, 0, rhs_begins.back() );
			
			job_stack.emplace_back( GEMM, a, b, c );
			
			for( auto rb_it=rhs_begins.begin(); cb_it != rhs_begins.end()-1; ++cb_it ){
				
				std::size_t rhs_begin = *rb_it;
				std::size_t rhs_end = *(rb_it+1);
				
				auto a = A.slice( diag_begin, diag_end, diag_begin, diag_end );
				auto b = B.slice( diag_begin, diag_end, rhs_begin, rhs_end );
				
				job_stack.emplace_back( TRSM, a, b );
			}
		}
	}
	else if( side == "R" && uplo == "L" ){
	
		std::vector<std::size_t> diag_begins = B.get_col_begins();
		std::vector<std::size_t> rhs_begins = B.get_row_begins();
		
		for( auto db_it=diag_begins.end()-1; db_it != diag_begins.begin(); --db_it ){
			
			std::size_t diag_begin = *(db_it-1);
			std::size_t diag_end = *db_it;
			
			auto a = B.slice( 0, rhs_begins.back(), diag_end, diag_begins.back() );
			auto b = A.slice( diag_end, diag_begins.back(), diag_begin, diag_end );
			auto c = B.slice( 0, rhs_begins.back(), diag_begin, diag_end );
			
			job_stack.emplace_back( GEMM, a, b, c );
			
			for( auto rb_it=rhs_begins.begin(); cb_it != rhs_begins.end()-1; ++cb_it ){
				
				std::size_t rhs_begin = *rb_it;
				std::size_t rhs_end = *(rb_it+1);
				
				auto a = A.slice( diag_begin, diag_end, diag_begin, diag_end );
				auto b = B.slice( rhs_begin, rhs_end, diag_begin, diag_end );
				
				job_stack.emplace_back( TRSM, a, b );
			}
		}
	}
	else {
		throw std::logic_error( "Invalid side and/or uplo" );
	}
	
	return job_stack;
}


void LUH_trsm( char side, char uplo, char diag, hMatrixInterface A, hMatrixInterface B ){

	std::list<TRSM_job_descriptor> job_stack;
	job_stack.emplace_back( TRSM, A, B );

	for( ; !job_stack.empty(); job_stack.pop_front() ){

		auto job = job_stack.front();

		if( job.type == GEMM ){
			// TODO wrap a (or b) in dense_mat - OR give iterator the logic to present as dense if pointing to of diag of LU block
			hMatrixGEMM( -1.0, job.a, job.b, job.c );
		}
		else /* job_type == TRSM */ {

			auto a = job.a;
			auto b = job.b;

			auto type_b = b.block_type();

			if( type_b == DenseMatrix ){
				Leaf_TRSM( side, uplo, diag, a, b );
			}
			else if ( type_b == LowRankMatrix ){
				// TODO wrap appropriate basis
				hMatrix b_wrapper;
				if( side == "L" ) b_wrapper.mat = b.left;
				else b_wrapper.mat = b.right;
				Leaf_TRSM( side, uplo, diag, a, b_wrapper );
			}
			else {
				std::list<TRSM_job_descriptor> new_jobs = queue_H_TRSM( side, uplo, a, b, s );
				// Push new jobs onto stack right behind 'front'
				auto insert_pos = job_stack.begin();
				++insert_pos;
				job_stack.insert( insert_pos, new_jobs.begin(), new_jobs.end() );
			}
		}
	}
}



void hMatrixTRSM( char side, char uplo, char diag, hMatrixInterface A, hMatrixInterface B ){

	std:vector<typename hMatrixInterface::TreeIterator> LU_blocks = A.get_LU_blocks();
	// Put blocks in the correct order
	std::sort( LU_blocks.begin(), LU_blocks.end(), 
		[]( const typename hMatrixInterface::TreeIterator& a, const typename hMatrixInterface::TreeIterator& b ){ 
			return a.row_begin < b.row_begin;
		}
	);
		
	if( side == "L" && uplo == "L" || side == "R" && uplo == "U" ){
		for( auto LU_it=LU_blocks.begin(); LU_it != LU_blocks.end(); ++LU_it ){
			auto diag_begin = LU_it->row_begin;
			auto diag_end = LU_it->row_end;
			
			auto lu = A.slice( diag_begin, diag_end, diag_begin, diag_end );
		
			if( side == "L" && uplo == "L" ){
				auto b = B.slice( diag_begin, diag_end, 0, B.ncols() );
				
				LUH_trsm( side, uplo, diag, lu, b );
				
				auto a = A.slice( diag_end, A.nrows(), diag_begin, diag_end );
				auto c = B.slice( diag_end, B.nrows(), 0, B.ncols() );
				hMatrixGEMM( -1.0, a, b, c );
			}
			else /*( side == "R" && uplo == "U" )*/ {
				auto a = B.slice( 0, B.nrows(), diag_begin, diag_end );
				
				LUH_trsm( side, uplo, diag, lu, a );
				
				auto b = A.slice( diag_begin, diag_end, diag_end, A.ncols() );
				auto c = B.slice( 0, B.nrows(), diag_end, B.ncols() );
				hMatrixGEMM( -1.0, a, b, c );
			}
		}
	}
	// Iterate in reverse order
	else if ( side == "L" && uplo == "U" || side == "R" && uplo == "L" ) {
		for( auto LU_it=LU_blocks.rbegin(); LU_it != LU_blocks.rend(); ++LU_it ){
			auto diag_begin = LU_it->row_begin;
			auto diag_end = LU_it->row_end;
			
			auto lu = A.slice( diag_begin, diag_end, diag_begin, diag_end );
		
			if( side == "L" && uplo == "U" ){
				auto a = A.slice( diag_begin, diag_end, diag_end, A.ncols() );
				auto b = B.slice( diag_end, B.nrows(), 0, B.ncols() );
				auto c = B.slice( diag_begin, diag_end, 0, B.ncols() );
				hMatrixGEMM( -1.0, a, b, c );
				
				LUH_trsm( side, uplo, diag, lu, c );
			}
			else /*( side == "R" && uplo == "L" )*/ {
				auto a = B.slice( 0, B.nrows(), diag_end, B.ncols() );
				auto b = A.slice( diag_end, A.nrows(), diag_begin, diag_end );
				auto c = B.slice( 0, B.nrows(), diag_begin, diag_end );
				hMatrixGEMM( -1.0, a, b, c );
				
				LUH_trsm( side, uplo, diag, lu, c );
			}
		}
	}
	else {
		throw std::logic_error( "Invalid side and/or uplo" );
	}
}


 

