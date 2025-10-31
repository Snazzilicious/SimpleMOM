
// TODO
// Hierarchy must not change in any routine
// subroutine error checking, debug and release, and early exiting
// Matrix multiply
//    All combinations
//        LR+LR should still check for same dataset
//    rSVD
///   Force truncation of insert into Low Rank blocks
// Distributed Memory
//    how to "send" slices
//    slices must stop at the in-core block level i.e. a loadable hMatrix is the leaf data of an OOC hMatrix
//        loading a slice needs to be aware that its not pointing to the whole matrix
//    Ref counting on disk?
// Distributed Disk
// API organization


void augment_low_rank( hMatrixInterface M, std::size_t new_rank ){
	// TODO
}

void LR_to_LR_GEMM( alpha, A, B, C ){
	hMatrix new_left_basis, new_right_basis, new_operand;
	
	if( A.block_type() == hMatrix::BlockType::LowRank ){
		new_left_basis.dense_mat = A.left;
		new_right_basis.dense_mat = zeros;
		new_operand = A.right;
		HHL_gemm( new_operand, B, new_right_basis ); // will not recall this function
	}
	else {
		new_left_basis.dense_mat = zeros;
		new_right_basis.dense_mat = B.right;
		new_operand = B.left;
		HHL_gemm( A, new_operand, new_left_basis ); // will not recall this function
	}
	
	// TODO copy new_left_basis.mat into c.left, unless is the same
	// TODO copy new_right_basis.mat into c.right, unless is the same
	rSVD( C );
}


void Leaf_GEMM( hMatrixInterface A, hMatrixInterface B, hMatrixInterface C ){
	// TODO
	auto type_a = a.block_type();
	auto type_b = b.block_type();
	auto type_c = c.block_type();
	
	if( type_c == DenseMatrix ){
		// going straight in
	
	}
	else {
		// augment C to receive result
	
		rSVD( c );
	}
}


template<hMatrixSlice>
struct GEMM_job_descriptor{
	hMatrixSlice a;
	hMatrixSlice b;
	hMatrixSlice c;
	GEMM_job_descriptor( a, b, c );
}


std::vector<std::size_t> get_partition_union( const std::vector<std::size_t>& part1, const std::vector<std::size_t>& part2 ){
	std::vector<std::size_t> part_union;
	std::set_union( part1.begin(), part1.end(), part2.begin(), part2.end(), std::back_inserter(part_union) );
	return part_union;
}

std::list<GEMM_job_descriptor> queue_H_GEMM( hMatrixInterface A, hMatrixInterface B, hMatrixInterface C ){
	// Get unions of partition points
	std::vector<std::size_t> row_begins = A.get_row_begins();
	std::vector<std::size_t> col_begins = B.get_col_begins();
	std::vector<std::size_t> inr_begins = get_partition_union( A.get_col_begins(), B.get_row_begins() );
	
	// Symbolically matmul all sub blocks
	std::list<GEMM_job_descriptor> job_stack;
	for( auto rb_it=row_begins.begin(); rb_it != row_begins.end()-1; ++rb_it ){
		for( auto cb_it=col_begins.begin(); cb_it != col_begins.end()-1; ++cb_it ){
			for( auto ib_it=inr_begins.begin(); ib_it != inr_begins.end()-1; ++ib_it ){
	
				std::size_t row_begin = *rb_it;
				std::size_t row_end = *(rb_it+1);
				std::size_t col_begin = *cb_it;
				std::size_t col_end = *(cb_it+1);
				std::size_t inr_begin = *ib_it;
				std::size_t inr_end = *(ib_it+1);
				
				auto a = A.slice( row_begin, row_end, inr_begin, inr_end );
				auto b = B.slice( inr_begin, inr_end, col_begin, col_end );
				auto c = C.slice( row_begin, row_end, col_begin, col_end );
			
				job_stack.emplace_back( a, b, c );
			}
		}
	}

	return job_stack;
}


void HHL_gemm( Scalar alpha, hMatrixInterface A, hMatrixInterface B, hMatrixInterface C ){
	
	std::list<GEMM_job_descriptor> job_stack;
	job_stack.emplace_back( A, B, C );
	
	for( ; !job_stack.empty(); job_stack.pop_front() ){
		auto job = job_stack.front();

		auto a = job.a;
		auto b = job.b;
		auto c = job.c;
		
		auto type_a = a.block_type();
		auto type_b = b.block_type();
		auto type_c = c.block_type();
		
		if( type_a != hMatrix::BlockType::Zero && type_b != hMatrix::BlockType::Zero ){
			continue;
		}
		// Special case of LowRank times non-LowRank assigned to LowRank or Zero
		else if( ( type_c == hMatrix::BlockType::Zero || type_c == hMatrix::BlockType::LowRank ) 
		     && (( type_a == hMatrix::BlockType::LowRank ) != ( type_b == hMatrix::BlockType::LowRank )) ){
			LR_to_LR_GEMM( alpha, a, b, c );
		}
		// Base case
		else if( type_a != hMatrix && type_b != hMatrix ){
			Leaf_GEMM( alpha, a, b, c );
		}
		// General case
		else {
			std::list<GEMM_job_descriptor> new_jobs = queue_H_GEMM( a, b, c );
			// Push new jobs onto stack right behind 'front'
			auto insert_pos = job_stack.begin();
			++insert_pos;
			job_stack.insert( insert_pos, new_jobs.begin(), new_jobs.end() );
		}
	}
}

template<Scalar,hMatrixInterface>
void hMatrixGEMM( Scalar alpha, hMatrixInterface A, hMatrixInterface B, hMatrixInterface C ){

	if ( !(A.nrows() == C.nrows() && B.ncols() == C.ncols() && A.ncols() == B.nrows()) )
		throw std::runtime_error("hMatrixGEMM: Incompatible operand dimensions.");

	if ( alpha == Scalar(0.0) ) return;
	
	std:vector<typename hMatrixInterface::TreeIterator> C_leaves = C.get_leaves();
	
	#pragma omp parallel for
	for( auto& C_leaf_it : C_leaves ){
		auto a = A.slice( C_leaf_it.row_begin, C_leaf_it.row_end, 0, A.ncols() );
		auto b = B.slice( 0, B.nrows(), C_leaf_it.col_begin, C_leaf_it.col_end );
		auto c = C.slice( C_leaf_it.row_begin, C_leaf_it.row_end, C_leaf_it.col_begin, C_leaf_it.col_end );
		HHL_gemm( alpha, a, b, c );
	}
}




template<Scalar,OOChMatrixInterface>
void hMatrixGEMM( Scalar alpha, OOChMatrixInterface A, OOChMatrixInterface B, OOChMatrixInterface C, std::size_t n_threads, MPI_comm Work_Comm ){

	// Get all leaves and associated rows/cols in C
	std:vector<typename OOChMatrixInterface::TreeIterator> C_icbs = C.get_in_core_blocks();
	
	// decide which procs do which leaves
	
	// For each icb
	//    get slices of A,B to multiply
	//    call gemm leaf c
	
	//load C
	for( std::size_t step=0; step < n_steps; ++step ){
		
	}
	// write C
	
}


template<Scalar,OOChMatrixInterface>
void hMatrixGEMM( Scalar alpha, OOChMatrixInterface A, OOChMatrixInterface B, OOChMatrixInterface C,
                  std::size_t n_threads, MPI_comm Work_Comm, std::unordered_map<std::string,std::vector<int>> proc_groups ){

	// TODO
}

