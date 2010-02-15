#ifndef OBJECTMATRIX_PEER_UTILS_HPP
#define OBJECTMATRIX_PEER_UTILS_HPP

#include <functional>
#include <string>

#include <Python.h>

#include <boost/python/object.hpp>
#include <boost/python/exception_translator.hpp>

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
        Tval_ const& get(T_ const& impl)
        {
            return (impl.*Vgetter_)();
        }

        void set(T_& impl, Tval_ const& v)
        {
            (impl.*Vsetter_)() = v;
        }
    };

    namespace detail {
        template<typename T1_, typename T2_, T1_(T2_::*Vfun_)()>
        struct range_from_range
        {
            static typename T1_::iterator begin(T2_& inst)
            {
                return (inst.*Vfun_)().begin();
            }

            static typename T1_::iterator end(T2_& inst)
            {
                return (inst.*Vfun_)().end();
            }
        };

        template<typename T1_, typename T2_, T1_(T2_::*Vfun_)() const>
        struct range_from_range<T1_, const T2_, Vfun_>
        {
            static typename T1_::const_iterator begin(T2_ const& inst)
            {
                return (inst.*Vfun_)().begin();
            }

            static typename T1_::const_iterator end(T2_ const& inst)
            {
                return (inst.*Vfun_)().end();
            }
        };
    } // namespace detail

    template<typename T1_, typename T2_, T1_(T2_::*Vfun_)() const>
    inline boost::python::object range_from_range()
    {
        return boost::python::range(
            &detail::range_from_range<T1_, const T2_, Vfun_>::begin,
            &detail::range_from_range<T1_, const T2_, Vfun_>::end);
    }

    template<typename T1_, typename T2_, T1_(T2_::*Vfun_)()>
    inline boost::python::object range_from_range()
    {
        return boost::python::range(
            &detail::range_from_range<T1_, T2_, Vfun_>::begin,
            &detail::range_from_range<T1_, T2_, Vfun_>::end);
    }
 
    template<typename T1_, typename T2_, T1_(T2_::*Vfun_)() const, typename Tpol_>
    inline boost::python::object range_from_range()
    {
        return boost::python::range<Tpol_>(
            &detail::range_from_range<T1_, const T2_, Vfun_>::begin,
            &detail::range_from_range<T1_, const T2_, Vfun_>::end);
    }

    template<typename T1_, typename T2_, T1_(T2_::*Vfun_)(), typename Tpol_>
    inline boost::python::object range_from_range()
    {
        return boost::python::range<Tpol_>(
            &detail::range_from_range<T1_, T2_, Vfun_>::begin,
            &detail::range_from_range<T1_, T2_, Vfun_>::end);
    }

} // namespace util

} // namespace peer

#endif /* OBJECTMATRIX_PEER_UTILS_HPP */
