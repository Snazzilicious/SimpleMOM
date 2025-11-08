
// TODO
// Leaf blocks - template on type
//    slice these
// Set item
//    matrix data insert vs hierarchy insert
//    need to identify and replace leaf nodes
//    LR should check if has same dataset, maybe
//    Hierarchy will not change in any routine
// types for OOC and distributed matrices need to all be mutually disjoint, or need to be different node types
// Constructing hMatrices
//    Assigning matrix data to Dense and Low rank blocks
//    concatenation
// come up with container for nodes
//    and replace node ptrs with iterators
// OOC variants
//    may want to generalize Slice to avoid duplication
//    this could be at python level
//    should load from arbitrary slice
// Need a generic tree traversal and node retrieval algorithm - see BGL
// abstract out index type
// const



void check_slice_limits( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end, std::size_t ncols, std::size_t nrows ){
	if( row_begin > row_end || row_end > nrows || col_begin > col_end || col_end > ncols )
		throw std::runtime_error("Invalid slice range.");
}


template<typename Scalar>
struct MatrixData {
	private:
		std::unique_ptr<Scalar> _data;
		int _begin_offset;
		int _nrows, _ncols;
		Layout _layout;
		int _ld;
	
	public:
		MatrixData( int rows, int cols ) : _nrows(nrows), _ncols(ncols), _begin_offset(0), _layout(COL_MAJOR), ld(nrows) {
			data = make_shared<Scalar>( _nrows*_ncols );
		}
		
		MatrixData( const MatrixData& m ) : _data(m._data) _nrows(m._nrows), _ncols(m._ncols), _begin_offset(m._begin_offset), _layout(m._layout), _ld(m._ld) {}
		
		enum class Layout { COL_MAJOR, ROW_MAJOR }
		
		Scalar* get_ptr(){ return _data.get() + _begin_offset; }
		int nrows(){ return _nrows; }
		int ncols(){ return _ncols; }
		Layout layout(){ return _layout; }
		int cblas_layout(){ return _layout; }
		int lapacke_layout(){ return _layout; }
		int ld(){ return _ld; }
		
		MatrixData slice( int row_begin, int row_end, int col_begin, int col_end ){
			check_slice_limits( row_begin, row_end, col_begin, col_end, _nrows, _ncols );
			
			MatrixData s(*this);
			
			s._begin_offset += layout == COL_MAJOR ? row_begin + col_begin*ld : row_begin*ld + col_begin ;
			s._nrows = row_end-row_begin;
			s._ncols = col_end-col_begin;
			
			return s;
		}
		
		MatrixData transpose(){
			MatrixData t(*this);
			
			t._layout = _layout == COL_MAJOR ? ROW_MAJOR : COL_MAJOR ;
			t._nrows = _ncols;
			t._ncols = _nrows;
		}
		
		Scalar& at( int i, int j ){
			check_slice_limits( i, i+1, j, j+1, _nrows, _ncols );
			
			return _layout == COL_MAJOR ? _data[ i + j*_ld ] : _data[ i*_ld + j ];
		}
};


struct Payload {
	private:
		BlockType _block_type;
	
	protected:
		Payload( BlockType type ) : _block_type(type) {}
	
	public:
		BlockType block_type(){ return _block_type; }
		virtual ~Payload(){}
		virtual std::size_t nrows()=0;
		virtual std::size_t ncols()=0;
};


struct ZeroBlock : public Payload {
	private:
		std::size_t _nrows,_ncols;
	
	public:
		ZeroBlock( std::size_t n_rows, std::size_t n_cols ) : Payload( BlockType::Zero ), _nrows(n_rows), _ncols(n_cols) {}
	
		std::size_t nrows() override { return _nrows; }
		std::size_t ncols() override { return _ncols; }
};


template<typename Scalar>
struct DenseBlock : public Payload {
	private:
		MatrixData<Scalar> mat;
	
	public:
		DenseBlock( std::size_t n_rows, std::size_t n_cols ) : mat(n_rows,n_cols) {}
		DenseBlock( const MatrixData& m ) : mat(m) {}
		
		MatrixData& mat(){ return mat; }
	
		std::size_t nrows() override { return mat.nrows(); }
		std::size_t ncols() override { return mat.ncols(); }
};


template<typename Scalar>
struct LowRankBlock : public Payload {
	private:
		MatrixData<Scalar> left,right;
	
	public:
		LowRankBlock( std::size_t n_rows, std::size_t n_cols, std::size_t rank ) : left(n_rows,rank), right(rank,n_cols) {}
		
		MatrixData& left(){ return left; }
		MatrixData& right(){ return right; }
		
		std::size_t nrows() override { return left.nrows(); }
		std::size_t ncols() override { return right.ncols(); }
		std::size_t rank(){ return left.ncols(); }
};


template<typename Scalar>
class ChildArray : public Payload {
	private:
		std::size_t _nrows,_ncols;
		MatrixData<hMatrix<Scalar>> blocks;
	
	public:
		ChildArray( std::size_t n_rows, std::size_t n_cols, std::size_t n_row_blocks, std::size_t n_col_blocks )
		 : _nrows(n_rows), _ncols(n_cols), blocks(n_row_blocks,n_col_blocks) {}
		 
		 hMatrix<Scalar>& at( std::size_t i, std::size_t j ){ return blocks.at(i,j); }
	
		std::size_t nrows() override { return _nrows; }
		std::size_t ncols() override { return _ncols; }
		
		std::size_t n_row_blocks() override { return blocks.nrows(); }
		std::size_t n_col_blocks() override { return blocks.ncols(); }
};




template<typename Scalar>
class hMatrix {
	private:
		std::unique_ptr<Payload> payload;
		
		void resize( std::size_t n_rows, std::size_t n_cols ){
			payload = std::make_unique<ZeroBlock>( n_rows, n_cols );
		}
	
	public:
		hMatrix( std::size_t n_rows, std::size_t n_cols ){
			payload = std::make_unique<ZeroBlock>( n_rows, n_cols );
		}
		
		enum class BlockType { Zero, Dense, LowRank, H };
		
		BlockType block_type(){ return payload->block_type(); }
		
		std::size_t nrows(){ return payload->nrows(); }
		std::size_t ncols(){ return payload->ncols(); }
		
		std::vector<std::size_t> get_row_begins(){
			std::vector<std::size_t> begins;
			if( this->block_type() == BlockType::H ){
				auto children = get_child_array();
				begins.resize( children.n_row_blocks()+1 );
				begins[0] = 0;
				for( std:size_t i=0; i<begins.size()-1; ++i )
					begins[i+1] = begins[i] + children->blocks.at( i, 0 ).nrows();
			}
			else {
				begins = { 0, nrows() };
			}
			return begins;
		}
		std::vector<std::size_t> get_col_begins(){
			std::vector<std::size_t> begins;
			if( this->block_type() == BlockType::H ){
				auto children = get_child_array();
				begins.resize( children.n_col_blocks()+1 );
				begins[0] = 0;
				for( std:size_t i=0; i<begins.size()-1; ++i )
					begins[i+1] = begins[i] + children->blocks.at( 0, i ).ncols();
			}
			else {
				begins = { 0, ncols() };
			}
			return begins;
		}
		
		
		
		DenseBlock<Scalar>& get_dense_block();
		LowRankBloc<Scalar>& get_low_rank_block();
		ChildArray<Scalar>& get_child_array();
		
		
		// find lowest node in the tree which completely contains the slice
		Slice slice( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end ){
			check_slice_limits( row_begin, row_end, col_begin, col_end, this->nrows(), this->ncols() );
			
			Node* node = this;
			bool fits_in_child = true;
			while( fits_in_child && node->block_type() == BlockType::H ){
				std::vector<std::size_t> row_begins = node->get_row_begins();
				std::vector<std::size_t> col_begins = node->get_col_begins();
				
				// find first child intersected by slice
				auto row_block_it = std::find_if( row_begins.begin(), row_begins.end(), 
					[row_begin](std::size_t block_begin){
						return block_begin<=row_begin;
					}
				);
				auto col_block_it = std::find_if( col_begins.begin(), col_begins.end(), 
					[col_begin](std::size_t block_begin){
						return block_begin<=col_begin;
					}
				);
				std::size_t row_block = row_block_it - row_begins.begin();
				std::size_t col_block = col_block_it - col_begins.begin();
				
				fits_in_child = row_end <= row_begins[row_block+1] && col_end <= col_begins[col_block+1] ;
				
				// slice child
				if( fits_in_child ){
					row_begin -= row_begins[row_block];
					row_end -= row_begins[row_block];
					col_begin -= col_begins[col_block];
					col_end -= col_begins[col_block];
					
					node = &(node->get_child_array().at(row_block, col_block ));
				}
			}
			
			return Slice( *node, row_begin, row_end, col_begin, col_end );
		}
		
		
		// for inserting hierarchical nodes
		void partition( const std::vector<std::size_t>& row_begins, const std::vector<std::size_t>& col_begins )
		{
			if( block_type() != BlockType::Zero )
				throw std::runtime_error( "Attempt to partition non-Zero block." ); // is possible to insert anything anywhere, just not worth implementing
			if( !std::is_sorted( row_begins.begin(), row_begins.end() ) || !std::is_sorted( col_begins.begin(), col_begins.end() ) )
				throw std::runtime_error( "Partition points are unsorted." );
			if( row_begins.back() > this->nrows() || col_begins.back() > this->ncols() )
				throw std::runtime_error( "Partition value(s) out of bounds." );
			
			// TODO account for 0 and end
			
			std::size_t n_rows = nrows();
			std::size_t n_cols = ncols();
			payload = std::make_unique<ChildArray>( n_rows, n_cols, row_begins.size()-1, col_begins.size()-1 );
			
			auto children = get_child_array();
			for( std::size_t i=0; i<row_begins.size()-1; ++i ){
				for( std::size_t j=0; j<col_begins.size()-1; ++j ){
					n_rows = row_begins[i+1]-row_begins[i];
					n_cols = col_begins[i+1]-col_begins[i];
					children.at( i,j ).resize( n_rows, n_cols );
				}
			}
		}
	
	
		// for inserting leaf data
		// Modifies tree structure, inserts exact values, returns Slice to output
		// always makes a copy of data to be inserted. To avoid this, insert empty nodes then assign to slices
		// assumes ownership of the matrix data pointed to by the 
		Slice insert( DenseBlock<Scalar> p ){
		Slice insert( LowRankBlock<Scalar> p ){
			// error checking
			
			// allocate node.payload as appropriate type
			
			// copy *p
		}
		
		// for assigning to leaves, preserving hierarchical structure
		// Preserves tree structure, may truncate
		void assign( Slice input, Slice output, Real tol=1e-5 );
};




class hMatrix::Slice {
	private:
		Node& root; // make this an iterator to whatever contains these - will require pointer / reference to tree object
		std::size_t rbegin, rend, cbegin, cend;
		
		
		static std::vector<std::size_t> begins_in_range( const std::vector<std::size_t>& begins, std::size_t begin, std::size_t end )
		{
			auto in_range_comp = [begin,end]( std::size_t v ){ return begin < v && v < end; };
			
			auto n_begins = std::count_if( begins.begin(), begins.end(), in_range_comp );
			
			n_begins += 2;
			std::vector<std::size_t> in_range_begins( n_begins, 0 );
			in_range_begins[ n_begins-1 ] = end;
			
			std::copy_if( begins.begin(), begins.end(), in_range_begins.begin(), in_range_comp );
		}
		
	public:
		Slice( Node& root_node, std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end )
			: root(root_node), rbegin(row_begin), rend(row_end), cbegin(col_begin), cend(col_end)
		{
			check_slice_limits( row_begin, row_end, col_begin, col_end, root_node->nrows(), root_node->ncols() );
		}
		
		BlockType block_type(){ return root->block_type(); }
		
		std::size_t nrows(){ return rend-rbegin; }
		std::size_t ncols(){ return cend-cbegin; }
		
		
		std::vector<std::size_t> get_row_begins(){
			return begins_in_range( root->get_row_begins(), rbegin, rend );
		}
		std::vector<std::size_t> get_col_begins(){
			return begins_in_range( root->get_col_begins(), cbegin, cend );
		}
		
		
		DenseBlockSlice get_dense_block();
		LowRankBlockSlice get_low_rank_block();
		
		
		// find lowest node in the tree which completely contains the slice
		Slice slice( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end ){
			check_slice_limits( row_begin, row_end, col_begin, col_end, this->nrows(), this->ncols() );
			
			row_begin += rbegin;
			row_end += rbegin;
			col_begin += cbegin;
			col_end += cbegin;
			
			return root.slice( row_begin, row_end, col_begin, col_end );
		}
};











// OOC Stuff



struct Payload {
	private:
		BlockType _block_type;
	
	public:
		Payload( BlockType type ) : _block_type(type) {}
		BlockType block_type(){ return _block_type; }
		virtual ~Payload(){}
};


template<typename Scalar>
struct OOCBlock : public Payload {
	FileDescriptor f;
};

class ChildArray : public Payload {
	std::vector<hMatrix::Node> blocks;
	int nrows, ncols;
	int layout;
	int ld();
};











