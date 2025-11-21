

#include<shared_ptr>
#include<limits>

template<typename Scalar>
struct MatrixData {
	private:
		std::shared_ptr<Scalar> _data;
		Scalar* _begin_ptr;
		std::size_t _nrows, _ncols;
		Layout _layout;
		std::size_t _ld;
	
	public:
		MatrixData();		
		MatrixData( std::size_t rows, std::size_t cols );
		MatrixData( const MatrixData& m );
		MatrixData( Scalar* data_begin, std::size_t n_rows, std::size_t n_cols, Layout layout, std::size_t stride ); // For borrowing data
		
		enum class Layout { COL_MAJOR, ROW_MAJOR }
		
		static int cblas_layout( Layout l );
		static int lapacke_layout( Layout l );
		template<typename Int_t>
		Int_t convert_to_int( std::size_t i ){
			if( i > static_cast<Index_t>(std::numeric_limits<Int_t>::max()) )
				throw std::runtime_error( "value too large to convert to requested type." );
			return std::static_cast<Index_t>(_nrows);
		}
		
		Scalar* get_ptr();
		
		std::size_t nrows();
		std::size_t ncols();
		Layout layout();
		std::size_t ld();
		
		MatrixData slice( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end );
		
		MatrixData transpose();
		
		Scalar& at( std::size_t i, std::size_t j );
};



template<typename Scalar>
MatrixData<Scalar>::MatrixData() : _nrows(0), _ncols(0), _begin_ptr(nullptr), _layout(Layout::COL_MAJOR), _ld(0), _data(nullptr) {}

template<typename Scalar>
MatrixData<Scalar>::MatrixData( std::size_t n_rows, std::size_t n_cols ) : _nrows(n_rows), _ncols(n_cols), _layout(Layout::COL_MAJOR), _ld(nrows) {
	data = make_shared<Scalar>( _nrows*_ncols );
	_begin_ptr = data.get();
}

template<typename Scalar>
MatrixData<Scalar>::MatrixData( const MatrixData& m ) 
: _data(m._data) _nrows(m._nrows), _ncols(m._ncols), _begin_ptr(m._begin_ptr), _layout(m._layout), _ld(m._ld) {}

template<typename Scalar>
MatrixData<Scalar>::MatrixData( Scalar* data_begin, std::size_t n_rows, std::size_t n_cols, Layout layout, std::size_t stride )
: _data(nullptr), _nrows(n_rows), _ncols(n_cols), _begin_ptr(data_begin), _layout(layout), _ld(stride) {} // For borrowing data

template<typename Scalar>
Scalar* MatrixData<Scalar>::get_ptr(){ return _begin_ptr; }

template<typename Scalar>
Index_t MatrixData<Scalar>::nrows(){ return _nrows; }
template<typename Scalar>
Index_t MatrixData<Scalar>::ncols(){ return _ncols; }

template<typename Scalar>
Layout MatrixData<Scalar>::layout(){ return _layout; }

template<typename Scalar>
Index_t MatrixData<Scalar>::ld(){ return _ld; }

template<typename Scalar>
MatrixData MatrixData<Scalar>::slice( std::size_t row_begin, std::size_t row_end, std::size_t col_begin, std::size_t col_end ){
	check_slice_limits( row_begin, row_end, col_begin, col_end, _nrows, _ncols );
	
	MatrixData s(*this);
	
	s._begin_ptr += layout == Layout::COL_MAJOR ? row_begin + col_begin*ld : row_begin*ld + col_begin ;
	s._nrows = row_end-row_begin;
	s._ncols = col_end-col_begin;
	
	return s;
}

template<typename Scalar>
MatrixData MatrixData<Scalar>::transpose(){
	MatrixData t(*this);
	
	t._layout = _layout == Layout::COL_MAJOR ? Layout::ROW_MAJOR : Layout::COL_MAJOR ;
	t._nrows = _ncols;
	t._ncols = _nrows;
}

template<typename Scalar>
Scalar& MatrixData<Scalar>::at( std::size_t i, std::size_t j ){
	check_slice_limits( i, i+1, j, j+1, _nrows, _ncols );
	
	return _layout == Layout::COL_MAJOR ? _data[ i + j*_ld ] : _data[ i*_ld + j ];
}




