#ifndef POSITION3_TRAITS_HPP
#define POSITION3_TRAITS_HPP

// This file containing the reference to the ecell4::Real3 and 
//   some template traits of ecell4::Real3.
//

#include <ostream>
#include <iomanip>
#include <functional>
#include <algorithm>

#include <ecell4/core/Real3.hpp>

#include <boost/array.hpp>
#include "utils/array_traits.hpp"
#include "linear_algebra.hpp"

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
struct is_vector<ecell4::Real3, N_>: public boost::mpl::true_ {};

template <>
struct element_type_of<ecell4::Real3>
{
    typedef ecell4::Real3::value_type type;
};

template<>
struct shape_position_type<ecell4::Real3>
{
    typedef ecell4::Real3 type;
};

template<>
struct shape_length_type<ecell4::Real3>
{
    typedef ecell4::Real3::value_type type;
};

inline ecell4::Real3 shape_position(ecell4::Real3 const &v)
{
    return v;
}

#endif /* POSITION3_TRAITS_HPP */
