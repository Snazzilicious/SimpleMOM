
#include "Ray.hpp"


RayList::iterator RayList::begin(){ return iterator( this ); }
RayList::iterator RayList::end(){ return begin()+size(); }
