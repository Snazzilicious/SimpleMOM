
// TODO
// OOC variants
//    may want to generalize Slice to avoid duplication
//    this could be at python level
//    should load from arbitrary slice
// Need a generic tree traversal and node retrieval algorithm - see BGL
// abstract out index type
// const
// this

// Assign / +=
	// d->d
	// d->lr : d -= lr, lr' = rSVD(d), lr += lr', lr = rSVD(lr), optionally d' = todense(lr) / d += lr, lr = rSVD(d) else d
	// lr->d
	// lr->lr
// rSVD
	// d, lr

// maybe Slice should be completely self sufficient, since it's supposed to act like a complete matrix itself
// hMatrix methods should only take a TreeIterator
// Yes I like this - see python

template<typename Scalar>
class hMatrix {
	private:
		std::size_t _min_block_size;
		
		std::vector<MatrixData<std::size_t>> connectivity; // connectivity, OR matrix indices
		
		std::vector<BlockType> types;
		
		std::vector<std::size_t> node_rows, node_cols;
		
		std::vector<MatrixData<Scalar>> matrices;
		
	public:
		enum class BlockType { Zero, Dense, LowRank, H };
		
		class TreeIterator;
		class Slice;
		
		// these work even for slices (assuming they align with boundaries)
		static std::size_t num_blocks( std::size_t min_block_size, std::size_t length );
		static std::size_t block_size( std::size_t min_block_size, std::size_t length, std::size_t block_index );
		static std::size_t block_index( std::size_t min_block_size, std::size_t length, std::size_t position );
		static std::size_t block_begin( std::size_t min_block_size, std::size_t length, std::size_t block_index );
		static bool on_block_boundary( std::size_t min_block_size, std::size_t length, std::size_t position );
		
		hMatrix( std::size_t n_rows, std::size_t n_cols, std::size_t min_block_size=256 );
		
		std::size_t nrows();
		std::size_t ncols();
		std::size_t min_block_size();
		
		TreeIterator root_node();
		Slice root_slice();
		
		void validate_iterator( const TreeIterator& node );
		
		std::size_t nrows( const TreeIterator& node );
		std::size_t ncols( const TreeIterator& node );
		BlockType block_type( const TreeIterator& node );
		
		// Main tree construction routines
		void partition( const TreeIterator& node, std::vector<std::size_t>& row_parts, std::vector<std::size_t>& col_parts );
		TreeIterator insert_dense( const TreeIterator& node );
		TreeIterator insert_lowrank( const TreeIterator& node );
		TreeIterator insert_dense( const TreeIterator& node, const MatrixData<Scalar>& D );
		TreeIterator insert_lowrank( const TreeIterator& node, const MatrixData<Scalar>& L, const MatrixData<Scalar>& R );
		
		// Main tree traversal routines
		std::size_t nrow_children( const TreeIterator& node );
		std::size_t ncol_children( const TreeIterator& node );
		TreeIterator get_child( const TreeIterator& node, std::size_t row_i, std::size_t col_i );
		
		std::vector<std::size_t> get_row_bounds( const TreeIterator& node );
		std::vector<std::size_t> get_col_bounds( const TreeIterator& node );
		static std::vector<std::size_t> bounds_in_range( const std::vector<std::size_t>& bounds, std::size_t begin, std::size_t end );
		
		// must be on block boundaries
		Slice slice( const TreeIterator& start_node, std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end );
		
		// TODO
		MatrixData<Scalar> get_dense_data( const TreeIterator& node );
		std::pair<MatrixData<Scalar>,MatrixData<Scalar>> get_lowrank_data( const TreeIterator& node );
		
		void serialize( char* buf ); // TODO
		
		hMatrix change_min_block_size( std::size_t new_min_block_size ); // TODO check if new size is compatible before changing ... maybe
}


class hMatrix::TreeIterator {
	private:
		hMatrix *parent_matrix;
		std::size_t _node_index;
		
		TreeIterator( hMatrix *parent_matrix, std::size_t node_ID );
	
	public:
		std::size_t nrows();
		std::size_t ncols();
		BlockType block_type();
		
		Slice root_slice();
		
		// Main tree construction routines
		void partition( const std::vector<std::size_t>& row_parts, const std::vector<std::size_t>& col_parts );
		void partition( std::vector<std::size_t>& row_parts, std::vector<std::size_t>& col_parts );
		TreeIterator insert_dense();
		TreeIterator insert_lowrank();
		TreeIterator insert_dense( const MatrixData<Scalar>& D );
		TreeIterator insert_lowrank( const MatrixData<Scalar>& L, const MatrixData<Scalar>& R );
		
		// Main tree traversal routines
		std::size_t nrow_children();
		std::size_t ncol_children();
		TreeIterator get_child( std::size_t row_i, std::size_t col_i );
		
		std::vector<std::size_t> get_row_bounds();
		std::vector<std::size_t> get_col_bounds();
		
		// must be on block boundaries
		Slice slice( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end );
}


class hMatrix::Slice {
	private:
		TreeIterator parent_node;
		std::size_t rbegin, rend, cbegin, cend;
		
		Slice( TreeIterator node, std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end );
		
	public:
		TreeIterator get_node();
		std::size_t nrows();
		std::size_t ncols();
		BlockType block_type();
		
		std::vector<std::size_t> get_row_bounds();
		std::vector<std::size_t> get_col_bounds();
		
		// must be on block boundaries
		Slice slice( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end );
		
		// TODO
		void plus_assign( const MatrixData<Scalar>& D );
		void plus_assign( const MatrixData<Scalar>& L, const MatrixData<Scalar>& R );
		
		// TODO
		MatrixData<Scalar> get_dense_data();
		std::pair<MatrixData<Scalar>,MatrixData<Scalar>> get_lowrank_data();
}


void check_slice_limits( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end, std::size_t ncols, std::size_t nrows ){
	if( row_begin > row_end || row_end > nrows || col_begin > col_end || col_end > ncols )
		throw std::runtime_error("Invalid slice range.");
}



static std::size_t hMatrix::num_blocks( std::size_t min_block_size, std::size_t length )
{
	return std::max( 1, length / min_block_size );
}

static std::size_t hMatrix::block_size( std::size_t min_block_size, std::size_t length, std::size_t block_index )
{
	std::size_t n_blocks = num_blocks( min_block_size, length );
	
	if( block_index >= n_blocks )
		throw std::runtime_error("Index of block for requested size exceeds num blocks.");
	
	if( length < min_block_size )
		return length;
	
	std::size_t remainder = length % min_block_size ;
	
	return min_block_size + std::static_cast<std::size_t>( block_index < remainder );
}

static bool hMatrix::block_index( std::size_t min_block_size, std::size_t length, std::size_t position )
{
	if( position >= length )
		throw std::runtime_error("Position exceeds length");
	
	std::size_t remainder = length % min_block_size ;
	
	std::size_t position_part1 = std::min( position, remainder * (min_block_size+1) );
	std::size_t position_part2 = position - position_part1;
	
	return ( position_part1 / (min_block_size+1) ) + ( position_part2 / min_block_size );
}

static std::size_t hMatrix::block_begin( std::size_t min_block_size, std::size_t length, std::size_t block_index )
{
	std::size_t n_blocks = num_blocks( min_block_size, length );
	
	if( block_index >= n_blocks )
		throw std::runtime_error("Index of block for requested size exceeds num blocks.");
		
	std::size_t remainder = length % min_block_size ;
	
	return min_block_size * block_index + std::min( remainder, block_index );
}

static bool hMatrix::on_block_boundary( std::size_t min_block_size, std::size_t length, std::size_t position )
{
	if( position > length )
		throw std::runtime_error("Position exceeds length");
	
	if( position == length )
		return true;
	
	std::size_t index = block_index( min_block_size, length, position );
	
	return position == block_begin( min_block_size, length, index );
}




hMatrix::hMatrix( std::size_t n_rows, std::size_t n_cols, std::size_t min_block_size ) : _min_block_size(min_block_size)
{	
	connectivity.emplace_back();
	types.push_back( BlockType::Zero );
	node_rows.push_back( n_rows );
	node_cols.push_back( n_cols );
}

std::size_t hMatrix::nrows(){ return nrows( root_node() ); }
std::size_t hMatrix::ncols(){ return ncols( root_node() ); }

std::size_t hMatrix::min_block_size(){ return _min_block_size; }

TreeIterator hMatrix::root_node(){ return TreeIterator( this, 0 ); }
Slice hMatrix::root_slice(){ return Slice( root_node(), 0, nrows(), 0, ncols() ); }.


void hMatrix::validate_iterator( const TreeIterator& node ){
	if( node._parent != this )
		throw std::runtime_error( "Iterator not associated with this hMatrix." );
	if( node._node_index >= connectivity.size() )
		throw std::runtime_error( "Iterator node index out of bounds." );
}

void hMatrix::validate_slice( const Slice& submat ){}


std::size_t hMatrix::nrows( const TreeIterator& node ){
	validate_iterator( node );
	return nrows[ node._node_index ];
}
std::size_t hMatrix::ncols( const TreeIterator& node ){
	validate_iterator( node );
	return ncols[ node._node_index ];
}
BlockType hMatrix::block_type( const TreeIterator& node ){
	validate_iterator( node );
	return types[ node._node_index ];
}


void hMatrix::partition( const TreeIterator& node, std::vector<std::size_t>& row_parts, std::vector<std::size_t>& col_parts )
{
	validate_iterator( node );
	
	std::size_t n_rows = this->nrows( node );
	std::size_t n_cols = this->ncols( node );
	std::size_t block_size = this->min_block_size();
	
	// ensure this is a zero block
	if( block_type( node ) != BlockType::Zero )
		throw std::runtime_error("Cannot partition non-zero block.");
	
	// ensure parts arrays are sorted and in bounds and on block boundaries
	if( !std::is_sorted( row_parts.begin(), row_parts.end() ) )
		throw std::runtime_error( "Row partition points are unsorted." );
	if( !std::is_sorted( col_parts.begin(), col_parts.end() ) )
		throw std::runtime_error( "Col partition points are unsorted." );
	if( row_begins.back() > n_rows || col_begins.back() > n_cols )
		throw std::runtime_error( "Partition value(s) out of bounds." );
	if( !all_of( row_parts.begin(), row_parts.end(), [n_rows,block_size](std::size_t i){ return on_block_boundary( block_size, n_rows, i ); } ) )
		throw std::runtime_error( "Row partition point(s) are unaligned with block boundaries." );
	if( !all_of( col_parts.begin(), col_parts.end(), [n_cols,block_size](std::size_t i){ return on_block_boundary( block_size, n_cols, i ); } ) )
		throw std::runtime_error( "Col partition point(s) are unaligned with block boundaries." );
	
	// add 0 and end if not present
	if( row_parts.front() != 0 )
		row_parts.insert( row_parts.begin(), 0 );
	if( row_parts.back() != n_rows )
		row_parts.push_back( n_rows );
	if( col_parts.front() != 0 )
		col_parts.insert( col_parts.begin(), 0 );
	if( col_parts.back() != n_cols )
		col_parts.push_back( n_cols );
	
	connectivity[ node._node_index ] = MatrixData<std::size_t>( row_parts.size()-1, col_parts.size()-1 );
	types[ node._node_index ] = BlockType::H;
	
	for( std::size_t i=0; i<row_parts.size()-1; ++i ){
		for( std::size_t j=0; j<col_parts.size()-1; ++j ){
			std::size_t n_rows = row_begins[i+1]-row_begins[i];
			std::size_t n_cols = col_begins[j+1]-col_begins[j];
			
			// create new zero block
			connectivity.emplace_back();
			types.push_back( BlockType::Zero );
			node_rows.push_back( n_rows );
			node_cols.push_back( n_cols );
			
			// add edge to new child block
			connectivity[ node._node_index ].at( i, j ) = connectivity.size()-1;
		}
	}
}


TreeIterator hMatrix::insert_dense( const TreeIterator& node ){
	validate_iterator( node );
	types[ node._node_index ] = BlockType::Dense;
}
TreeIterator hMatrix::insert_lowrank( const TreeIterator& node ){
	validate_iterator( node );
	types[ node._node_index ] = BlockType::LowRank;
}

TreeIterator hMatrix::insert_dense( const TreeIterator& node, const MatrixData<Scalar>& D ){
	insert_dense( node );
	matrices.emplace_back( D ); // TODO check dimensions
	connectivity[ node._node_index ] = MatrixData<std::size_t>( 1,1 );
	connectivity[ node._node_index ].at( 0,0 ) = matrices.size()-1 ;
}
TreeIterator hMatrix::insert_lowrank( const TreeIterator& node, const MatrixData<Scalar>& L, const MatrixData<Scalar>& R ){
	insert_lowrank( node );
	matrices.emplace_back( L ); // TODO check these dimensions
	matrices.emplace_back( R );
	connectivity[ node._node_index ] = MatrixData<std::size_t>( 1,2 );
	connectivity[ node._node_index ].at( 0,0 ) = matrices.size()-2 ;
	connectivity[ node._node_index ].at( 0,1 ) = matrices.size()-1 ;
}



		
std::size_t nrow_children( const TreeIterator& node ){
	validate_iterator( node );
	return connectivity[ node._node_index ].nrows();
}
std::size_t ncol_children( const TreeIterator& node )
	validate_iterator( node );
	return connectivity[ node._node_index ].ncols();
		
TreeIterator hMatrix::get_child( const TreeIterator& node, std::size_t i, std::size_t j ){
	validate_iterator( node );
	
	if( block_type( node ) != BlockType::H )
		throw std::runtime_error( "Cannot get child of non-H node." );
	if( i >= nrow_children( node ) || j >= ncol_children( node ) )
		throw std::runtime_error( "Child indices out of bounds" );
	
	return TreeIterator( this, connectivity[node._node_index].at(i,j) );
}


		
std::vector<std::size_t> hMatrix::get_row_bounds( const TreeIterator& node ){
	validate_iterator( node );
	
	std::vector<std::size_t> bounds;
	
	if( block_type( node ) == BlockType::H )
	{
		std::size_t n_children = nrow_children( node );
		bounds.resize( n_children + 1 );
		bounds[0] = 0;
		for( std:size_t i=0; i<n_children; ++i ){
			TreeIterator child = get_child( node, i, 0 );
			bounds[i+1] = bounds[i] + nrows( child );
		}
	}
	else 
	{
		bounds = { 0, nrows( node ) };
	}
	return bounds;
}

std::vector<std::size_t> hMatrix::get_col_bounds( const TreeIterator& node ){
	validate_iterator( node );
	
	std::vector<std::size_t> bounds;
	
	if( block_type( node ) == BlockType::H )
	{
		std::size_t n_children = ncol_children( node );
		bounds.resize( n_children + 1 );
		bounds[0] = 0;
		for( std:size_t i=0; i<n_children; ++i ){
			TreeIterator child = get_child( node, 0, i );
			bounds[i+1] = bounds[i] + ncols( child );
		}
	}
	else 
	{
		bounds = { 0, ncols( node ) };
	}
	return bounds;
}



// find lowest node in the tree which completely contains the slice
Slice slice( const TreeIterator& start_node, std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end ){
	validate_iterator( start_node );
	
	check_slice_limits( row_begin, row_end, col_begin, col_end, this->nrows(), this->ncols() );
	
	// TODO ensure on block boundaries
	
	TreeIterator node = start_node;
	bool fits_in_child = true;
	while( fits_in_child && block_type( node ) == BlockType::H ){
		std::vector<std::size_t> row_bounds = get_row_bounds( node );
		std::vector<std::size_t> col_bounds = get_col_bounds( node );
		
		// find first child intersected by slice
		auto row_block_it = std::find_if( row_bounds.begin(), row_bounds.end(), 
			[row_begin](std::size_t block_begin){
				return block_begin<=row_begin;
			}
		);
		auto col_block_it = std::find_if( col_bounds.begin(), col_bounds.end(), 
			[col_begin](std::size_t block_begin){
				return block_begin<=col_begin;
			}
		);
		std::size_t row_block = row_block_it - row_bounds.begin();
		std::size_t col_block = col_block_it - col_bounds.begin();
		
		fits_in_child = row_end <= row_bounds[row_block+1] && col_end <= col_bounds[col_block+1] ;
		
		// slice child
		if( fits_in_child ){
			row_begin -= row_bounds[row_block];
			row_end -= row_bounds[row_block];
			col_begin -= col_bounds[col_block];
			col_end -= col_bounds[col_block];
			
			node = get_child( node, row_block, col_block );
		}
	}
	
	return Slice( node, row_begin, row_end, col_begin, col_end );
}

//------------------------------------------------------------------------------------------------------------------------------

hMatrix::TreeIterator::TreeIterator( hMatrix *parent_matrix, std::size_t node_ID ) : _parent(parent_matrix), _node_index(node_ID) {}

// These all just call the parent matrix's routine with *this passed in

std::size_t hMatrix::TreeIterator::nrows(){ return parent_matrix->nrows( *this ); }
std::size_t hMatrix::TreeIterator::ncols(){ return parent_matrix->ncols( *this ); }
BlockType hMatrix::TreeIterator::block_type(){ return parent_matrix->block_type( *this ); }

Slice hMatrix::TreeIterator::root_slice(){ return slice( 0, nrows(), 0, ncols() ); }

void partition( std::vector<std::size_t>& row_parts, std::vector<std::size_t>& col_parts )
{
	parent_matrix->partition( *this, row_parts, col_parts );
}
TreeIterator insert_dense(){ return parent_matrix->insert_dense( *this ); }
TreeIterator insert_lowrank(){ return parent_matrix->insert_lowrank( *this ); }
TreeIterator insert_dense( const MatrixData<Scalar>& D ){ return parent_matrix->insert_dense( *this, D ); }
TreeIterator insert_lowrank( const MatrixData<Scalar>& L, const MatrixData<Scalar>& R ){ return parent_matrix->insert_lowrank( *this, L, R ); }

std::size_t nrow_children(){ return parent_matrix->nrow_children( *this ); }
std::size_t ncol_children(){ return parent_matrix->ncol_children( *this ); }
TreeIterator get_child( std::size_t row_i, std::size_t col_i ){ return parent_matrix->get_child( *this, i, j ); }

std::vector<std::size_t> get_row_bounds(){ return parent_matrix->get_row_bounds( *this ); }
std::vector<std::size_t> get_col_bounds(){ return parent_matrix->get_col_bounds( *this ); }

Slice slice( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end )
{
	return parent_matrix->slice( *this, row_begin, row_end, col_begin, col_end );
}

//------------------------------------------------------------------------------------------------------------------------------


hMatrix::Slice::Slice( TreeIterator node, std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end )
	: parent_node(node), rbegin(row_begin), rend(row_end), cbegin(col_begin), cend(col_end) {}

TreeIterator get_node(){ return parent_node; }

// These mostly just call the parent matrix's routine with *this passed in

std::size_t Slice::nrows(){ return rend-rbegin; }
std::size_t Slice::ncols(){ return cend-cbegin; }
BlockType Slice::block_type(){ return parent_node.block_type(); }


std::vector<std::size_t> bounds_in_range( const std::vector<std::size_t>& bounds, std::size_t begin, std::size_t end )
{
	auto in_range_comp = [begin,end]( std::size_t v ){ return begin < v && v < end; };
	
	auto n_bounds = std::count_if( bounds.begin(), bounds.end(), in_range_comp );
	
	n_bounds += 2;
	std::vector<std::size_t> in_range_bounds( n_bounds, 0 );
	in_range_bounds[ n_bounds-1 ] = end;
	
	std::copy_if( bounds.begin(), bounds.end(), in_range_bounds.begin(), in_range_comp );
	
	return in_range_bounds;
}


std::vector<std::size_t> Slice::get_row_bounds(){
	auto row_bounds = parent_node.get_row_bounds();
	return bounds_in_range( row_bounds, rbegin, rend );
}

std::vector<std::size_t> Slice::get_col_bounds()
	auto col_bounds = parent_node.get_col_bounds();
	return bounds_in_range( col_bounds, cbegin, cend );
}

Slice Slice::slice( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end ){
	return parent_node.slice( row_begin+rbegin, row_end+rbegin, col_begin+cbegin, col_end+cbegin );
}

void plus_assign( const MatrixData<Scalar>& D ){ return parent_node.parent_matrix->plus_assign( *this, D ); }
void plus_assign( const MatrixData<Scalar>& L, const MatrixData<Scalar>& R ){ return parent_node.parent_matrix->plus_assign( *this, L, R ); }













