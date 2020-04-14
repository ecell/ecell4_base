#ifndef POSITION3_TRAITS_HPP
#define POSITION3_TRAITS_HPP

// This file containing the reference to the ecell4::Real3 and 
//   some template traits of ecell4::Real3.
//

#include <ostream>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <array>

#include <ecell4/core/Real3.hpp>

#include "utils/array_traits.hpp"
#include "linear_algebra.hpp"

namespace ecell4
{
namespace egfrd
{
template<typename T_>
struct shape_position_type
{
    struct argument_is_not_a_shape;
    static const std::size_t x = sizeof(argument_is_not_a_shape);
};

template<typename T_>
struct shape_length_type
{
    typedef typename element_type_of<typename shape_position_type<typename boost::remove_cv<T_>::type >::type>::type type;
};

template <std::size_t N_>
struct is_vector<Real3, N_>: public std::true_type {};

template <>
struct element_type_of<Real3>
{
    typedef Real3::value_type type;
};

template<>
struct shape_position_type<Real3>
{
    typedef Real3 type;
};

template<>
struct shape_length_type<Real3>
{
    typedef Real3::value_type type;
};

inline ::ecell4::Real3 shape_position(::ecell4::Real3 const &v)
{
    return v;
}

} // egfrd
} // ecell4
#endif /* POSITION3_TRAITS_HPP */
