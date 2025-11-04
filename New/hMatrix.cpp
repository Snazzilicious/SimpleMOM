
// TODO
// Leaf blocks - template on type
//    poylmorphic shared_ptrs
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


class hMatrix {
	private:
		class Node;
		class Slice;
		std::shared_ptr<Node> root; // make this an iterator to whatever contains these
		
		// the graph
		
		// hash_table of DenseBlocks
		
		// hash_table of LowRankBlocks
		
		// hash_table of children configurations
	
	public:
		hMatrix( std::size_t n_rows, std::size_t n_cols ){
			root = std::make_shared<Node>( n_rows, n_cols );
		}
		
		
		enum class BlockType { Zero, Dense, LowRank, H };
		
		Slice slice( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end ){
			Slice s( root, 0, root->nrows(), 0, root->ncols() );
			return s.slice( row_begin, row_end, col_begin, col_end );
		}
		Slice slice(){
			return slice( 0, root->nrows(), 0, root->ncols() );
		}
		
		
		void partition( Slice position, const std::vector<std::size_t>& row_begins, const std::vector<std::size_t>& col_begins ){
			if( position.block_type() != BlockType::Zero )
				throw std::runtime_error( "Can only insert into Zero block." ); // technically can insert anything anywhere, just not worth implementing
			
			
		}
		
		
		// Modifies tree structure, inserts exact values, returns Slice to output
		// always makes a copy of data to be inserted. To avoid this, insert empty nodes then assign to slices
		Slice insert( Slice position, payload_t p ){
			// partition first
			
			// 
		}
		
		// Preserves tree structure, may truncate
		void assign( Slice input, Slice output, Real tol=1e-5 );
};


struct Payload {
	private:
		BlockType _block_type;
	
	public:
		Payload( BlockType type ) : _block_type(type) {}
		BlockType block_type(){ return _block_type; }
		virtual ~Payload(){}
};

template<typename Scalar>
struct MatrixData {
	std::shared_ptr<Scalar> data;
	int nrows, ncols;
	int layout;
	int ld();
};

template<typename Scalar>
struct DenseBlock : public Payload {
	MatrixData<Scalar> D;
};

template<typename Scalar>
struct LowRankBlock : public Payload {
	MatrixData<Scalar> U,V;
};

class ChildArray : public Payload {
	std::vector<hMatrix::Node> blocks;
	int nrows, ncols;
	int layout;
	int ld();
};

template<typename Scalar>
class hMatrix::Node {
	private:
		std::size_t _nrows,_ncols;
		std::shared_ptr<Payload> payload;
	
	public:
		Node( std::size_t n_rows, std::size_t n_cols ) : _nrows(n_rows), _ncols(n_cols) {
			payload = std::make_shared<Payload>( BlockType::Zero );
		}
		
		BlockType block_type(){ return payload->block_type(); }
		
		std::size_t nrows(){ return _nrows; }
		std::size_t ncols(){ return _ncols; }
		
		std::vector<std::size_t> get_row_begins(){
			std::vector<std::size_t> begins;
			if( this->block_type() == BlockType::H ){
				auto children = get_child_array();
				begins.resize( children->nrows+1 );
				begins[0] = 0;
				std::size_t ld = children->layout == ROW_MAJOR ? 1 : children->ld();
				for( std:size_t i=0; i<begins.size()-1; ++i )
					begins[i+1] = begins[i] + children->blocks[ i*ld ].nrows();
			}
			else {
				begins.resize(2);
				begins[0] = 0;
				begins[1] = _nrows;
			}
			return begins;
		}
		std::vector<std::size_t> get_col_begins(){
			std::vector<std::size_t> begins;
			if( this->block_type() == BlockType::H ){
				auto children = get_child_array();
				begins.resize( children->ncols+1 );
				begins[0] = 0;
				std::size_t ld = children->layout == COL_MAJOR ? 1 : children->ld();
				for( std:size_t i=0; i<begins.size()-1; ++i )
					begins[i+1] = begins[i] + children->blocks[ i*ld ].ncols();
			}
			else {
				begins.resize(2);
				begins[0] = 0;
				begins[1] = _ncols;
			}
			return begins;
		}
		
		std::shared_ptr<DenseBlock> get_dense_block();
		std::shared_ptr<LowRankBlock> get_low_rank_block();
		std::shared_ptr<ChildArray> get_child_array();
};



template<typename Scalar>
struct MatrixDataSlice {
	std::shared_ptr<Scalar> data;
	Scalar *ptr;
	int nrows, ncols;
	int layout;
	int ld;
};

template<typename Scalar>
struct DenseBlockSlice {
	MatrixDataSlice<Scalar> D;
};

template<typename Scalar>
struct LowRankBlockSlice {
	MatrixDataSlice<Scalar> U,V;
};



void check_slice_limits( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end std::size_t ncols, std::size_t nrows ){
	if( row_begin > row_end || row_end > nrows || col_begin > col_end || col_end > ncols )
		throw std::runtime_error("Invalid slice range.");
}


class hMatrix::Slice {
	private:
		std::size_t rbegin, rend, cbegin, cend;
		std::shared_ptr<Node> root; // make this an iterator to whatever contains these - will require pointer / reference to tree object
	
	public:
		Slice( std::shared_ptr<Node> root_node, std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end )
			: root(root_node), rbegin(row_begin), rend(row_end), cbegin(col_begin), cend(col_end)
		{
			check_slice_limits( row_begin, row_end, col_begin, col_end, root_node->nrows(), root_node->ncols() );
		}
		
		BlockType block_type(){ return root->block_type(); }
		
		std::size_t nrows(){ return rend-rbegin; }
		std::size_t ncols(){ return cend-cbegin; }
		
		std::vector<std::size_t> get_row_begins(){
			std::vector<std::size_t> root_begins = root->get_row_begins();
			std::vector<std::size_t> begins = {0};
			for( auto v : root_begins )
				if( rbegin < v && v < rend )
					begins.push_back( v-rbegin );
			begins.push_back( rend-rbegin );
			return begins;
		}
		std::vector<std::size_t> get_col_begins(){
			std::vector<std::size_t> root_begins = root->get_col_begins();
			std::vector<std::size_t> begins = {0};
			for( auto v : root_begins )
				if( cbegin < v && v < cend )
					begins.push_back( v-cbegin );
			begins.push_back( cend-cbegin );
			return begins;
		}
		
		
		DenseBlockSlice get_dense_block();
		LowRankBlockSlice get_low_rank_block();
		
		
		// find lowest node in the tree which completely contains the slice
		Slice slice( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end ){
			check_slice_limits( row_begin, row_end, col_begin, col_end, this->nrows(), this->ncols() );
			
			row_begin += rbegin;
			row_end += rend;
			col_begin += cbegin;
			col_end += cend;
			
			auto node = root;
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
					
					node = node->children[row_block][col_block];
				}
			}
			
			return Slice( node, row_begin, row_end, col_begin, col_end );
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











