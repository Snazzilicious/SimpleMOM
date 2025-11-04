
// TODO
// Tri Solves
// MPI OOC
//    need matrix concatenation for solving across OOC part boundaries
//    Should probably be recursive to minimize IO
// Dist Disk


void Leaf_TRSM( char side, char uplo, IntIt ipiv, hMatrixInterface A, hMatrixInterface B ){

	if ( A.block_type() != hMatrix::BlockType::Dense || B.block_type() != hMatrix::BlockType::Dense )
		throw std::runtime_error("hMatrixTRSM: Invalid operand types.");
	
	// Get pointers of dense data, (and strides and layouts)
	auto A_ptr = A.data();
	int A_m = A.nrows();
	int A_n = A.ncols();
	int lda = A.ld();
	int A_layout = A.layout;
	
	auto B_ptr = B.data();
	int B_m = B.nrows();
	int B_n = B.ncols();
	int ldb = B.ld();
	int B_layout = B.layout;
	
	int *p = ipiv.data();
	
	Scalar alpha = 1.0;
	
	if( side == 'L' && uplo == 'L' ){
		// B = P^-1 B
		int err = LAPACKE_?laswp( B_layout, B_n, B_ptr, ldb, 1, B_m, ipiv, 1 );
		if( err ) throw std::runtime_error( "LASWP failed in Leaf_TRSM" );
		
		// B = L^-1 B
		A_uplo = A_layout == B_layout ? CblasLower : CblasUpper;
		A_trans = A_layout == B_layout ? CblasNoTrans : CblasTrans;
		cblas_trsm( B_layout, CblasLeft, A_uplo, A_trans, CblasUnit, A_m, B_n, &alpha, A_ptr, lda, B_ptr, ldb );
	}
	else if( side == 'R' && uplo == 'L' ){
		// B = B L^-1
		A_uplo = A_layout == B_layout ? CblasLower : CblasUpper;
		A_trans = A_layout == B_layout ? CblasNoTrans : CblasTrans;
		cblas_trsm( B_layout, CblasRight, A_uplo, A_trans, CblasUnit, A_m, B_n, &alpha, A_ptr, lda, B_ptr, ldb );
		
		// B = B P^-1
		order = B_layout == COL_MAJOR ? LAPACK_ROW_MAJOR : LAPACK_COL_MAJOR;
		int err = LAPACKE_?laswp( order, B_m, B_ptr, ldb, 1, B_n, ipiv, -1 );
		if( err ) throw std::runtime_error( "LASWP failed in Leaf_TRSM" );
	}
	else if( side == 'L' && uplo == 'U' ){
		// B = U^-1 B
		A_uplo = A_layout == B_layout ? CblasUpper : CblasLower;
		A_trans = A_layout == B_layout ? CblasNoTrans : CblasTrans;
		cblas_trsm( B_layout, CblasLeft, A_uplo, A_trans, CblasNonUnit, A_m, B_n, &alpha, A_ptr, lda, B_ptr, ldb );
	}
	else /* side == 'R' && uplo == 'U' */ {
		// B = B U^-1
		A_uplo = A_layout == B_layout ? CblasUpper : CblasLower;
		A_trans = A_layout == B_layout ? CblasNoTrans : CblasTrans;
		cblas_trsm( B_layout, CblasRight, A_uplo, A_trans, CblasNonUnit, A_m, B_n, &alpha, A_ptr, lda, B_ptr, ldb );
	}
}


void LUH_trsm( char side, char uplo, IntIt ipiv, hMatrixInterface A, hMatrixInterface B ){
	// Get B leaves
	std::vector<typename hMatrixInterface::TreeIterator> B_leaves = B.get_leaves();
	
	if( side == "L" ){
		// get column partitions
		std::set<std::size_t> unique_col_begins = {B.ncols()};
		for( const auto& it : B_leaves )
			unique_col_begins.insert( it.col_begin );
		std::vector<std::size_t> col_begins( unique_col_begins.size() );
		std::copy( unique_col_begins.begin(), unique_col_begins.end(), col_begins.begin() );
		
		// Solve on each partition
		for( auto cb_it=col_begins.begin(); cb_it != col_begins.end()-1; ++cb_it ){

			std::size_t col_begin = *cb_it;
			std::size_t col_end = *(cb_it+1);
			
			auto b = B.slice( 0, B.nrows(), col_begin, col_end );
			
			if( b.block_type() == hMatrix::BlockType::LowRank ){
				// TODO
				hMatrix b_wrapper;
				b_wrapper.dense = b.left ;
				Leaf_TRSM( side, uplo, diag, A, b_wrapper );
			}
			else if( b.block_type() == hMatrix::BlockType::Dense ){
				Leaf_TRSM( side, uplo, diag, A, b );
			}
			else if( b.block_type() == hMatrix::BlockType::H ){
				hMatrix dense_b = b.todense();
				Leaf_TRSM( side, uplo, diag, A, dense_b );
				b.assign( dense_b.slice() );
			}
			// continue if b is zero type
		}
	}
	else {
		// get row partitions
		std::set<std::size_t> unique_row_begins = {B.nrows()};
		for( const auto& it : B_leaves )
			unique_row_begins.insert( it.row_begin );
		std::vector<std::size_t> row_begins( unique_row_begins.size() );
		std::copy( unique_row_begins.begin(), unique_row_begins.end(), row_begins.begin() );
		
		// Solve on each partition
		for( auto rb_it=row_begins.begin(); rb_it != row_begins.end()-1; ++rb_it ){
	
			std::size_t row_begin = *rb_it;
			std::size_t row_end = *(rb_it+1);
			
			auto b = B.slice( row_begin, row_end, 0, B.ncols() );
			
			if( b.block_type() == hMatrix::BlockType::LowRank ){
				// TODO
				hMatrix b_wrapper;
				b_wrapper.dense = b.right
				Leaf_TRSM( side, uplo, diag, A, b_wrapper );
			}
			else if( b.block_type() == hMatrix::BlockType::Dense ){
				Leaf_TRSM( side, uplo, diag, A, b );
			}
			else if( b.block_type() == hMatrix::BlockType::H ){
				hMatrix dense_b = b.todense();
				Leaf_TRSM( side, uplo, diag, A, dense_b );
				b.assign( dense_b.slice() );
			}
			// continue if b is zero type
		}
	}
}


struct TRSM_job_descriptor {
	hMatrixInterface a,b,c;
	IntIt p;
};


std::list<TRSM_job_descriptor> queue_H_TRSM( char side, char uplo, IntIt ipiv, hMatrixInterface A, hMatrixInterface B ){
	
	std::vector<std::size_t> diag_begins = A.row_begins();
	
	std::list<TRSM_job_descriptor> job_stack;
	if( side == "L" && uplo == "L" || side == "R" && uplo == "U" ){
		for( auto db_it=diag_begins.begin(); db_it != diag_begins.end()-1; ++db_it ){
			
			std::size_t diag_begin = *db_it;
			std::size_t diag_end = *(db_it+1);
			
			auto lu = A.slice( diag_begin, diag_end, diag_begin, diag_end );
			auto p = ipiv + diag_begin;
			
			if( side == "L" && uplo == "L" ){
				auto b = B.slice( diag_begin, diag_end, 0, B.ncols() );
				
				job_stack.emplace_back( TRSM, p, lu, b );
				
				auto a = A.slice( diag_end, A.nrows(), diag_begin, diag_end );
				auto c = B.slice( diag_end, B.nrows(), 0, B.ncols() );
				
				job_stack.emplace_back( GEMM, a, b, c );
			}
			else /*( side == "R" && uplo == "U" )*/ {
				auto a = B.slice( 0, B.nrows(), diag_begin, diag_end );
				
				job_stack.emplace_back( TRSM, p, lu, a );
				
				auto b = A.slice( diag_begin, diag_end, diag_end, A.ncols() );
				auto c = B.slice( 0, B.nrows(), diag_end, B.ncols() );
				
				job_stack.emplace_back( GEMM, a, b, c );
			}
		}
	}
	// Iterate in reverse order
	else if ( side == "L" && uplo == "U" || side == "R" && uplo == "L" ) {
		for( auto db_it=diag_begins.end()-1; db_it != diag_begins.begin(); --db_it ){
			
			std::size_t diag_begin = *(db_it-1);
			std::size_t diag_end = *db_it;

			auto lu = A.slice( diag_begin, diag_end, diag_begin, diag_end );
			auto p = ipiv + diag_begin;
		
			if( side == "L" && uplo == "U" ){
				auto a = A.slice( diag_begin, diag_end, diag_end, A.ncols() );
				auto b = B.slice( diag_end, B.nrows(), 0, B.ncols() );
				auto c = B.slice( diag_begin, diag_end, 0, B.ncols() );
				
				job_stack.emplace_back( GEMM, a, b, c );
				job_stack.emplace_back( TRSM, p, lu, c );
			}
			else /*( side == "R" && uplo == "L" )*/ {
				auto a = B.slice( 0, B.nrows(), diag_end, B.ncols() );
				auto b = A.slice( diag_end, A.nrows(), diag_begin, diag_end );
				auto c = B.slice( 0, B.nrows(), diag_begin, diag_end );

				job_stack.emplace_back( GEMM, a, b, c );
				job_stack.emplace_back( TRSM, p, lu, c );
			}
		}
	}
	else {
		throw std::logic_error( "Invalid side and/or uplo" );
	}
	
	return job_stack;
}


template<IntIt,hMatrixInterface>
void hMatrixTRSM( char side, char uplo, IntIt ipiv, hMatrixInterface A, hMatrixInterface B ){

	std::list<TRSM_job_descriptor> job_stack;
	job_stack.emplace_back( TRSM, ipiv, A, B );
	
	for( ; !job_stack.empty(); job_stack.pop_front() ){

		auto job = job_stack.front();
		
		auto p = job.p;
		auto a = job.a;
		auto b = job.b;
		auto c = job.c;
		
		if( job.type == GEMM ){
			hMatrixGEMM( -1.0, a, b, c );
		}
		else /* job.type == TRSM */ {
			// Base case
			if( a.block_type() == hMatrix::BlockType::Dense ){
				LUH_trsm( side, uplo, diag, p, a, b );
			}
			// Invalid case
			else if( a.block_type() == hMatrix::BlockType::Zero || a.block_type() == hMatrix::BlockType::LowRank ){
				throw std::runtime_error("Rank-deficient diagonal block in TRSM");
			}
			// General case
			else {
				std::list<TRSM_job_descriptor> new_jobs = queue_H_TRSM( side, uplo, p, a, b );
				// Push new jobs onto stack right behind 'front'
				auto insert_pos = job_stack.begin();
				++insert_pos;
				job_stack.insert( insert_pos, new_jobs.begin(), new_jobs.end() );
			}
		}
	}
} 

