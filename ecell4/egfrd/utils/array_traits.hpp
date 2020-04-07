#ifndef ECELL4_EGFRD_UTILS_ARRAY_TRAITS_HPP
#define ECELL4_EGFRD_UTILS_ARRAY_TRAITS_HPP

namespace ecell4
{
namespace egfrd
{
template< typename T_ >
struct element_type_of{};

// template< typename T_, std::size_t N_ >
// struct element_type_of< T_[N_] >
// {
//     typedef T_ type;
// };
//
// template< typename T_, std::size_t N_ >
// struct element_type_of< std::array< T_, N_ > >
// {
//     typedef T_ type;
// };
//
// template< typename T_, typename Talloc_ >
// struct element_type_of< boost::multi_array< T_, 1, Talloc_ > >
// {
//     typedef T_ type;
// };
} // ecell4
} // egfrd
#endif// ECELL4_EGFRD_UTILS_ARRAY_TRAITS_HPP
