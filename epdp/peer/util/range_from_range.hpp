#ifndef PEER_UTIL_RANGE_FROM_RANGE_HPP
#define PEER_UTIL_RANGE_FROM_RANGE_HPP

#include <functional>
#include <boost/python/object.hpp>
#include <boost/python/exception_translator.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/mpl/if.hpp>
#include "utils/range.hpp"

#include "peer/wrappers/iterator/stl_iterator_wrapper.hpp"
#include "peer/wrappers/range/stl_container_wrapper.hpp"

namespace peer { namespace util {

namespace detail {
    template<typename T1_, typename T2_, T1_(T2_::*Vfun_)(), typename Trcg_, bool B1_>
    struct range_from_range_impl
    {
        typedef typename boost::remove_reference<T1_>::type range_type;
        typedef typename boost::mpl::if_<
            boost::is_const<T2_>,
            typename boost::range_const_iterator<range_type>::type,
            typename boost::range_iterator<range_type>::type>::type
                result_type;

        typedef peer::wrappers::stl_iterator_wrapper<result_type, boost::python::object, Trcg_> wrapper_type;

        static PyObject* create_wrapper(boost::python::back_reference<T2_&> backref)
        {
            wrapper_type::__class_init__(typeid(result_type).name());
            return wrapper_type::create((backref.get().*Vfun_)(),
                    backref.source());
        };

        static boost::python::object create()
        {
            return boost::python::make_function(&create_wrapper);
        }
    };

    template<typename T1_, typename T2_, T1_(T2_::*Vfun_)(), typename Trcg_>
    struct range_from_range_impl<T1_, T2_, Vfun_, Trcg_, true>
    {
        typedef typename boost::remove_reference<T1_>::type range_type;
        typedef typename boost::mpl::if_<
            boost::is_const<T2_>,
            typename boost::range_const_iterator<range_type>::type,
            typename boost::range_iterator<range_type>::type>::type
                result_type;

        struct holder
        {
            T1_ operator*() const
            {
                return range_;
            };

            holder(T1_ range, boost::python::object py)
                : range_(range), py_(py) {}

            T1_ range_;
            boost::python::object py_;
        };

        typedef peer::wrappers::stl_container_wrapper<
                range_type, holder,
                typename boost::mpl::if_<
                    boost::is_const<T2_>,
                    peer::wrappers::default_policy_generator<peer::wrappers::default_immutable_container_wrapper_policy>,
                    peer::wrappers::default_policy_generator<peer::wrappers::default_container_wrapper_policy> >::type,
                Trcg_> wrapper_type;

        static PyObject* create_wrapper(boost::python::back_reference<T2_&> backref)
        {
            wrapper_type::__class_init__(typeid(result_type).name());
            return wrapper_type::create(
                    holder(
                        (backref.get().*Vfun_)(),
                        backref.source()));
        };

        static boost::python::object create()
        {
            return boost::python::make_function(&create_wrapper);
        }
    };

    template<typename T1_, typename T2_, T1_(T2_::*Vfun_)(), typename Trcg_ = boost::python::return_by_value>
    struct range_from_range: range_from_range_impl<T1_, T2_, Vfun_, Trcg_, check_range_iterator_category<typename boost::remove_reference<T1_>::type, boost::random_access_traversal_tag>::value> {};

} // namespace detail

template<typename T1_, typename T2_, T1_(T2_::*Vfun_)() const>
inline boost::python::object range_from_range()
{
    return detail::range_from_range<T1_, const T2_, Vfun_>::create();
}

template<typename T1_, typename T2_, T1_(T2_::*Vfun_)()>
inline boost::python::object range_from_range()
{
    return detail::range_from_range<T1_, T2_, Vfun_>::create();
}

template<typename T1_, typename T2_, T1_(T2_::*Vfun_)() const, typename Tpol_>
inline boost::python::object range_from_range()
{
    return detail::range_from_range<T1_, const T2_, Vfun_, Tpol_>::create();
}

template<typename T1_, typename T2_, T1_(T2_::*Vfun_)(), typename Tpol_>
inline boost::python::object range_from_range()
{
    return detail::range_from_range<T1_, T2_, Vfun_, Tpol_>::create();
}

} } // namespace peer::util

#endif /* PEER_UTIL_RANGE_FROM_RANGE_HPP */
