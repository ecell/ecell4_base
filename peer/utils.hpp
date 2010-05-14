#ifndef OBJECTMATRIX_PEER_UTILS_HPP
#define OBJECTMATRIX_PEER_UTILS_HPP

#include <functional>
#include <string>

#include <Python.h>

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

namespace peer {

namespace util
{
    inline boost::python::object pass_through(boost::python::object o)
    {
        return o;
    }

    template<typename T_>
    inline PyObject* pystr_from_repr(const T_* obj)
    {
        std::ostringstream s;
        s << (*obj);
        std::string res(s.str());
        return PyString_FromStringAndSize(res.data(), res.size());
    }

    namespace functor
    {
        template<typename T_>
        struct to_python_converter_fun
            : public std::unary_function<T_, boost::python::object>
        {
            typedef T_ argument_type;
            typedef boost::python::object result_type;

            result_type operator()(const argument_type& src)
            {
                return boost::python::object(src);
            }
        };

        template<typename Talloc_>
        class default_initializer_fun
            : public std::unary_function<typename Talloc_::reference, void>
        {
        public:
            typedef typename Talloc_::reference argument_type;
            typedef void result_type;

        public:
            default_initializer_fun(Talloc_& alloc): alloc_(alloc) {}

            void operator()(argument_type ptr)
            {
                new(alloc_.address(ptr)) typename Talloc_::value_type();
            }

        private:
            Talloc_& alloc_;
        };

        template<typename Talloc_>
        class destructor_fun
            : public std::unary_function<typename Talloc_::reference, void>
        {
        public:
            typedef typename Talloc_::reference argument_type;
            typedef void result_type;

        public:
            destructor_fun(Talloc_& alloc): alloc_(alloc) {}

            void operator()(argument_type ptr)
            {
                typedef typename Talloc_::value_type value_type;
                ptr.~value_type();
            }

        private:
            Talloc_& alloc_;
        };
    } // namespace functor

    template<typename Tnative_, typename Tconverter_>
    inline void to_native_converter()
    {
        boost::python::converter::registry::push_back(
                &Tconverter_::convertible,
                reinterpret_cast<
                        boost::python::converter::constructor_function>(
                            &Tconverter_::construct),
                boost::python::type_id<Tnative_>());
    }

    template<typename Tnative_, typename Tconverter_>
    inline void to_native_lvalue_converter()
    {
        boost::python::converter::registry::insert(
                &Tconverter_::convert,
                boost::python::type_id<Tnative_>(),
                &Tconverter_::expected_pytype);
    }

    static void std_exception_translator( const std::exception& exc )
    {
      PyErr_SetString( PyExc_RuntimeError, exc.what() );
    }

    inline void register_std_exception_translator()
    {
        boost::python::register_exception_translator< std::exception >(
            std_exception_translator );
    }

    template<typename T_, typename Tval_,
            Tval_ const&(T_::*Vgetter_)() const,
            Tval_ &(T_::*Vsetter_)()>
    struct reference_accessor_wrapper
    {
        static Tval_ const& get(T_ const& impl)
        {
            return (impl.*Vgetter_)();
        }

        static void set(T_& impl, Tval_ const& v)
        {
            (impl.*Vsetter_)() = v;
        }
    };

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

} // namespace util

} // namespace peer

#endif /* OBJECTMATRIX_PEER_UTILS_HPP */
