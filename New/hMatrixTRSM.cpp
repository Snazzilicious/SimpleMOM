
// TODO
// Tri Solves
// MPI OOC
//    need matrix concatenation for solving across OOC part boundaries
// Dist Disk


void Leaf_TRSM( char side, char uplo, char diag, hMatrixInterface A, hMatrixInterface B ){

	if ( !(A.block_type() == LU_Block && B.block_type() == DenseMatrix) )
		throw std::runtime_error("hMatrixTRSM: Invalid operand types.");
	
	// Get pointers of dense data, (and strides and layouts)
	
	if( side == 'L' && uplo == 'L' ){
		LAPACKE_?laswp( matrix_layout, n, a, lda, k1, k2, ipiv, incx );
		trsm( 'Left', 'Lower', 'No transpose', 'Unit', n, nrhs, one, a, lda, b, ldb );
	}
	else if( side == 'R' && uplo == 'L' ){
		LAPACKE_?laswp( (matrix_layout+1) % 2, n, a, lda, k1, k2, ipiv, -1 );
		trsm( 'Right', 'Lower', 'Transpose', 'Unit', n, nrhs, one, a, lda, b, ldb );
	}
	else if( side == 'L' && uplo == 'U' ){
		trsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n, nrhs, one, a, lda, b, ldb );
	}
	else {
		trsm( 'Right', 'Upper', 'No transpose', 'Non-unit', n, nrhs, one, a, lda, b, ldb );
	}
}


void LUH_trsm( char side, char uplo, char diag, hMatrixInterface A, hMatrixInterface B ){
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


 

