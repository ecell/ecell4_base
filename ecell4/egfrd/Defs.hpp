#if !defined( EGFRD_DEFS_HPP )
#define EGFRD_DEFS_HPP

#include <cstddef>
#include <ecell4/core/types.hpp>

namespace ecell4
{
namespace egfrd
{

// typedef double Real;
// typedef long int Integer;
typedef ecell4::Real Real;
typedef ecell4::Integer Integer;
typedef unsigned long int UnsignedInteger;
typedef size_t Index;

// stringifiers.  see preprocessor manual
#define XSTR( S ) STR( S )
#define STR( S ) #S

#define IGNORE_RETURN (void)

} // egfrd
} // ecell4
#endif // EGFRD_DEFS_HPP
