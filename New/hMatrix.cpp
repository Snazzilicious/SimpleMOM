
// TODO
// Leaf blocks - template on type
// Set item
//    matrix data insert vs hierarchy insert
//    need to identify and replace leaf nodes
//    LR should check if has same dataset, maybe
//    Hierarchy will not change in any routine
// How to know when pointing to leaf node? extra hMatrix node? add "type" labels? - leaning towards both
//    types for OOC and distributed matrices need to all be mutually disjoint, or need to be different node types
// Constructing hMatrices
//    Assigning matrix data to Dense and Low rank blocks
//    Shared pointers to nodes + shared_from_this
//    How to create on the stack?
// come up with container for nodes
//    and replace node ptrs with iterators
// OOC variants
//    may want to generalize Slice to avoid duplication
// Need a generic tree traversal and node retrieval algorithm - see BGL


class hMatrix {
	private:
		class Node;
		class Slice;
		std::shared_ptr<Node> root; // make this an iterator to whatever contains these
	
	public:
		enum class BlockType { Zero, Dense, LowRank, H };
		
		Slice slice( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end ){
			Slice whole_matrix( root, 0, root->nrows(), 0, root->ncols() );
			return whole_matrix.slice( row_begin, row_end, col_begin, col_end );
		}
		Slice slice(){
			return slice( 0, root->nrows(), 0, root->ncols() );
		}
};

template<typename Scalar>
class hMatrix::Node {
	private:
		std::size_t _nrows,_ncols;
		std::shared_ptr<Node> children;
		std::shared_ptr<MatrixData> dense;
		std::shared_ptr<MatrixData> left;
		std::shared_ptr<MatrixData> right;
		BlockType _block_type;
	
	public:
		
		BlockType block_type(){ return _block_type; }
		
		std::size_t nrows(){ return _nrows; }
		std::size_t ncols(){ return _ncols; }
		
		std::vector<std::size_t> get_row_begins(){
			std::vector<std::size_t> begins = {0};
			for( col_child )
				begins.push_back( begins.back()+col_child.nrows() );
			return begins;
		}
		std::vector<std::size_t> get_col_begins(){
			std::vector<std::size_t> begins = {0};
			for( row_child )
				begins.push_back( begins.back()+row_child.ncols() );
			return begins;
		}
};

void check_slice_limits( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end std::size_t ncols, std::size_t nrows ){
	if( row_begin > row_end || row_end > nrows || col_begin > col_end || col_end > ncols )
		throw std::runtime_error("Invalid slice range.");
}

class hMatrix::Slice {
	private:
		std::size_t rbegin, rend, cbegin, cend;
		std::shared_ptr<Node> root; // make this an iterator to whatever contains these
	
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
		
		// find lowest node in the tree which completely contains the slice
		Slice slice( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end ){
			check_slice_limits( row_begin, row_end, col_begin, col_end, this->nrows(), this->ncols() );
			
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
