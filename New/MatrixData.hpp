



template<typename Scalar>
struct MatrixData {
	private:
		std::shared_ptr<Scalar> _data;
		Scalar* _begin_offset;
		int _nrows, _ncols;
		Layout _layout;
		int _ld;
	
	public:
		MatrixData() : _nrows(0), _ncols(0), _begin_offset(0), _layout(COL_MAJOR), _ld(0), _data(nullptr) {}
		
		MatrixData( int rows, int cols ) : _nrows(nrows), _ncols(ncols), _begin_offset(0), _layout(COL_MAJOR), _ld(nrows) {
			data = make_shared<Scalar>( _nrows*_ncols );
		}
		
		MatrixData( const MatrixData& m ) 
		: _data(m._data) _nrows(m._nrows), _ncols(m._ncols), _begin_offset(m._begin_offset), _layout(m._layout), _ld(m._ld) {}
		
		MatrixData( Scalar* data_begin, int n_rows, int n_cols, Layout layout, int stride ) : _data(nullptr) {} // For borrowing data
		
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
