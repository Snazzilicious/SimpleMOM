
// TODO
// MPI OOC
// Dist Disk

void Leaf_GETRF( IntIt piv, hMatrixInterface A ){
	auto A_ptr = A.data;
	int m = A.nrows;
	int n = A.ncols;
	int lda = A.ld();
	int matrix_layout = A.layout;
	
	int *p = ipiv.data();
	
	if( m != n )
		throw std::runtime_error( "Non-square matrix in Leaf_GETRF: ("+ std::to_string(m) + ", " + std::to_string(n) + ")" );
	
	int err = LAPACKE_getrf( matrix_layout, m, m, A_ptr, lda, p );
	
	if( err )
		throw std::runtime_error( "Error in Leaf_GETRF: " + to_string(err) );
}


struct GETRF_job_descriptor {
	enum JobType { GETRF, TRSM_GEMM } type;
	IntIt p;
	hMatrixInterface a,b,c,d;
	GETRF_job_descriptor( IntIt piv, hMatrixInterface A );
	GETRF_job_descriptor( IntIt piv, hMatrixInterface A, hMatrixInterface B, hMatrixInterface C, hMatrixInterface D );
}


template<IntIt,hMatrixInterface>
void hMatrixGETRF( IntIt piv, hMatrixInterface A ){

	std::list<GETRF_job_descriptor> job_stack;
	job_stack.emplace_back( GETRF, piv, A );
	
	for( ; !job_stack.empty(); job_stack.pop_front() ){

		auto job = job_stack.front();
		
		auto p = job.p;
		auto a = job.a;
		auto b = job.b;
		auto c = job.c;
		auto d = job.d;
		
		if( job.type == TRSM_GEMM ){
			hMatrixTRSM( "L", "L", p, a, b );
			hMatrixTRSM( "R", "U", p, a, c );
			hMatrixGEMM( -1.0, c, b, d );
		}
		else /* job.type == GETRF */ {
			// Base case
			if( a.block_type() == hMatrix::BlockType::Dense ){
				Leaf_GETRF( p, a );
			}
			// Invalid case
			else if( a.block_type() == hMatrix::BlockType::Zero || a.block_type() == hMatrix::BlockType::LowRank ){
				throw std::runtime_error("Rank-deficient diagonal block in GETRF");
			}
			// General case
			else {
				std::vector<std::size_t> diag_begins = a.row_begins();
				
				for( auto db_it=diag_begins.begin(); db_it != diag_begins.end()-1; ++db_it ){
					
					std::size_t diag_begin = *db_it;
					std::size_t diag_end = *(db_it+1);
					
					auto new_p = p + diag_begin;
					auto new_a = a.slice( diag_begin, diag_end, diag_begin, diag_end );
					auto new_b = a.slice( diag_begin, diag_end, diag_end, a.ncols() );
					auto new_c = a.slice( diag_end, a.nrows(), diag_begin, diag_end );
					auto new_d = a.slice( diag_end, a.nrows(), diag_end, a.ncols() );
					
					GETRF_job_descriptor new_getrf( GETRF, new_p, new_a );
					GETRF_job_descriptor new_trsm_gemm( TRSM_GEMM, new_p, new_a, new_b, new_c, new_d );
					
					// Push new jobs onto stack right behind 'front'
					auto insert_pos = job_stack.begin();
					++insert_pos;
					insert_pos = job_stack.insert( insert_pos, new_trsm_gemm );
					insert_pos = job_stack.insert( insert_pos, new_getrf );
				}
			}
		}
	}
}
