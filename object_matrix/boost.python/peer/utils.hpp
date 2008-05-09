#ifndef OBJECTMATRIX_PEER_UTILS_HPP
#define OBJECTMATRIX_PEER_UTILS_HPP

#include <boost/python/object.hpp>
#include <functional>
#include <string>
#include <Python.h>

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
    void to_native_converter()
    {
        boost::python::converter::registry::push_back(
                &Tconverter_::convertible,
                reinterpret_cast<
                        boost::python::converter::constructor_function>(
                            &Tconverter_::construct),
                boost::python::type_id<Tnative_>());
    }
} // namespace util

} // namespace peer

#endif /* OBJECTMATRIX_PEER_UTILS_HPP */
