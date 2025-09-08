
#ifndef ITERATORS_HPP
#define ITERATORS_HPP

#include<iterator>

// For accessing rows (columns) of a 2D data set stored in row-major (column-major) format.
// Is much like a StrideIterator, but the dereference operator returns an iterator
// to the row (column) instead of the element pointed to.
template<typename Iterator>
class Iterator2D {
	private:
		Iterator pos;
		std::size_t stride;
	
	public:
		using difference_type = typename std::iterator_traits<Iterator>::difference_type ;
		using value_type = Iterator ;
		using pointer = Iterator ;
		using reference = Iterator ;
		using iterator_category = std::random_access_iterator_tag ;
		
		Iterator2D( Iterator begin, std::size_t step_size ) : pos(begin), stride(step_size) {}
		
		Iterator2D operator++(){ pos+=stride; return *this; }
		Iterator2D operator+(difference_type n){ Iterator2D tmp = *this; tmp.pos += stride*n; return tmp; }
		Iterator2D& operator+=(difference_type n){ pos += stride*n; return *this; }
		friend difference_type operator-(const Iterator2D& a, const Iterator2D& b){ return (a.pos - b.pos)/a.stride; }
		
		friend bool operator==( const Iterator2D& a, const Iterator2D& b ){ return a.pos == b.pos; }
		friend bool operator!=( const Iterator2D& a, const Iterator2D& b ){ return !(a == b); }
		
		pointer operator[](difference_type n){ return *(operator+(n)); }
		pointer operator*(){ return pos; }
		pointer operator->(){ return pos; }
};

// For accessing elements at constant intervals in an array
template<typename Iterator>
class StrideIterator {
	private:
		Iterator pos;
		std::size_t stride;
	
	public:
		using difference_type = typename std::iterator_traits<Iterator>::difference_type ;
		using value_type = typename std::iterator_traits<Iterator>::value_type ;
		using pointer = Iterator ;
		using reference = typename std::iterator_traits<Iterator>::reference ;
		using iterator_category = std::random_access_iterator_tag ;
		
		StrideIterator( Iterator begin, std::size_t step_size ) : pos(begin), stride(step_size) {}
		
		StrideIterator operator++(){ pos+=stride; return *this; }
		StrideIterator operator+(difference_type n){ StrideIterator tmp = *this; tmp.pos += stride*n; return tmp; }
		StrideIterator& operator+=(difference_type n){ pos += stride*n; return *this; }
		friend difference_type operator-(const StrideIterator& a, const StrideIterator& b){ return (a.pos - b.pos)/a.stride; }
		
		friend bool operator==( const StrideIterator& a, const StrideIterator& b ){ return a.pos == b.pos; }
		friend bool operator!=( const StrideIterator& a, const StrideIterator& b ){ return !(a == b); }
		
		reference operator[](difference_type n){ return *(operator+(n)); }
		reference operator*(){ return *pos; }
		pointer operator->(){ return pos; }
};

// For partitioning arrays or other sets into multiple groups e.g. for balanced parallel processing.
// This allocates equal, contiguous partitions to each group.
template<typename Iterator>
void distribute_evenly( Iterator begin, Iterator end, std::size_t n_groups, std::size_t group_id, Iterator& group_begin, Iterator& group_end ){
	std::size_t n_per_group = (end-begin) / n_groups ;
	std::size_t remainder = (end-begin) % n_groups ;
	
	group_begin = begin + group_id*n_per_group + std::min( group_id, remainder );
	group_end = group_begin + n_per_group + (int)( group_id < remainder );
}

// For partitioning arrays or other sets into multiple groups e.g. for balanced parallel processing.
// This allocates according to a round-robin scheme.
template<typename Iterator>
std::pair<StrideIterator<Iterator>,StrideIterator<Iterator>> distribute_round_robin( Iterator begin, Iterator end, std::size_t n_groups, std::size_t group_id ){
	Iterator first = begin + group_id;
	std::size_t n_this_group = ( (end-first) + (n_groups-1) ) / n_groups ;
	
	return std::make_pair( StrideIterator( first, n_groups ), StrideIterator( first, n_groups )+n_this_group );
}



// Interface to strided data that makes it act like a normal, iterable container such as an array or vector.
template<typename Iterator>
class StridedView {
	private:
		StrideIterator<Iterator> pos;
		std::size_t _size;

	public:
		StridedView( Iterator begin, std::size_t step_size, std::size_t length=0 ) : pos(begin,step_size), _size(length) {}
		
		typename std::iterator_traits<Iterator>::reference operator[](typename std::iterator_traits<Iterator>::difference_type n){ return *(pos+n); }
		
		StrideIterator<Iterator> begin(){ return pos; }
		StrideIterator<Iterator> end(){ return pos+_size; }
		
		std::size_t size(){ return _size; }
		
		template<typename Iterable>
		void operator=(const Iterable& a){
			std::copy( a.begin(),a.end(), pos );
		}
		
		template<typename Iterable>
		operator Iterable(){
			Iterable a(_size);
			std::copy( begin(), end(), a.begin() );
			return a;
		}
};

// Like an Iterator2D, but dereferencing returns a StridedView instead of an iterator.
// Is useful for iterating over the rows (columns) of column-major (row-major) data.
template<typename Iterator>
class ViewStridedIterator {
	private:
		Iterator pos;
		std::size_t stride;
		std::size_t strideView_size;
	
	public:
		using difference_type = typename std::iterator_traits<Iterator>::difference_type ;
		using value_type = StridedView<Iterator> ;
		using pointer = Iterator ;
		using reference = StridedView<Iterator> ;
		using iterator_category = std::random_access_iterator_tag ;
		
		ViewStridedIterator( Iterator begin, std::size_t step_size, std::size_t strideView_length=0 ) : pos(begin), stride(step_size), strideView_size(strideView_length) {}
		
		ViewStridedIterator operator++(){ ++pos; return *this; }
		ViewStridedIterator operator+(difference_type n) const { ViewStridedIterator tmp = *this; tmp.pos += n; return tmp; }
		ViewStridedIterator& operator+=(difference_type n){ pos += n; return *this; }
		friend difference_type operator-(const ViewStridedIterator& a, const ViewStridedIterator& b){ return a.pos - b.pos; }
		
		friend bool operator==( const ViewStridedIterator& a, const ViewStridedIterator& b ){ return a.pos == b.pos; }
		friend bool operator!=( const ViewStridedIterator& a, const ViewStridedIterator& b ){ return !(a == b); }
		
		reference operator[](difference_type n) const { return *(operator+(n)); }
		reference operator*() const { return StridedView<Iterator>( pos, stride, strideView_size ); }
		pointer operator->() const { return pos; }
};




template<typename T>
class CountingIterator {
	private:
		T pos;
	
	public:
		using difference_type = T ;
		using value_type = T ;
		using pointer = T ;
		using reference = T ;
		using iterator_category = std::random_access_iterator_tag ;
		
		CountingIterator( T begin=0 ) : pos(begin) {}
		
		CountingIterator operator++(){ ++pos; return *this; }
		CountingIterator operator+(difference_type n){ CountingIterator tmp = *this; tmp.pos += n; return tmp; }
		CountingIterator& operator+=(difference_type n){ pos += n; return *this; }
		friend difference_type operator-(const CountingIterator& a, const CountingIterator& b){ return a.pos - b.pos; }
		
		friend bool operator==( const CountingIterator& a, const CountingIterator& b ){ return a.pos == b.pos; }
		friend bool operator!=( const CountingIterator& a, const CountingIterator& b ){ return !(a == b); }
		
		reference operator[](difference_type n){ return *(operator+(n)); }
		reference operator*(){ return pos; }
		pointer operator->(){ return pos; }
};




#include<unordered_map>
#include<numeric>

template<typename Iterator, typename T=typename std::iterator_traits<Iterator>::value_type>
std::unordered_map<T,std::pair<Iterator,Iterator>> get_value_ranges( Iterator begin, Iterator end ){ // range is assumed sorted
	
	std::unordered_map<T,std::size_t> counts;
	counts = std::transform_reduce( begin, end, counts,
		[]( std::unordered_map<T,std::size_t> a, std::unordered_map<T,std::size_t> b )
		{
			for( auto val_cnt : a )
				b[val_cnt->first] += val_cnt->second;
			return b;
		},
		[](T v){ return std::unordered_map<T,std::size_t>{{v,1}}; } );
	
	std::unordered_map<T,std::pair<Iterator,Iterator>> ranges;
	for( auto pos=begin; pos != end; pos += counts[*pos] )
		ranges[*pos] = std::make_pair( pos, pos + counts[*pos] );
	
	return ranges;
}

template<typename Iterator, typename UnaryOp, typename in_t=typename std::iterator_traits<Iterator>::value_type, typename out_t=typename std::invoke_result<UnaryOp,in_t>>
std::unordered_map<out_t,std::pair<Iterator,Iterator>> get_transformed_value_ranges( Iterator begin, Iterator end, UnaryOp f ){ // range is assumed sorted
	
	std::unordered_map<out_t,std::size_t> counts;
	counts = std::transform_reduce( begin, end, counts,
		[]( std::unordered_map<out_t,std::size_t> a, std::unordered_map<out_t,std::size_t> b )
		{
			for( auto val_cnt : a )
				b[val_cnt->first] += val_cnt->second;
			return b;
		},
		[&f](in_t v){ return std::unordered_map<out_t,std::size_t>{{f(v),1}}; } );
	
	std::unordered_map<out_t,std::pair<Iterator,Iterator>> ranges;
	for( auto pos=begin; pos != end; pos += counts[f(*pos)] )
		ranges[*pos] = std::make_pair( pos, pos + counts[f(*pos)] );
	
	return ranges;
}

template<typename T, typename Iterator>
void get_range_iterators( Iterator& begin, Iterator& end, std::unordered_map<T,std::pair<Iterator,Iterator>> ranges, T val ){
	auto val_begin_end = ranges.find( val );
	if( val_begin_end != ranges.end() ){
		begin = val_begin_end->second.first;
		end = val_begin_end->second.second;
	}
}


#endif /* ITERATORS_HPP */


