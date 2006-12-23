#if !defined( __DEFS_HPP )
#define __DEFS_HPP

#include <vector>
#include <boost/multi_array.hpp>


typedef double Real;
typedef int Int;
typedef size_t Index;


typedef std::vector< Real > RealVector;
//typedef boost::multi_array< Real, 1, boost::pool_allocator<Real> > 
//RealArray;
typedef boost::multi_array<Real, 2>
Real2DArray;
typedef boost::multi_array<Real, 3>
Real3DArray;
typedef boost::multi_array<Real, 4>
Real4DArray;


#endif // __DEFS_HPP
