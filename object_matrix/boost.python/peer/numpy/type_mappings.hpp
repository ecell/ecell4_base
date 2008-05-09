#ifndef OBJECTMATRIX_PEER_NUMPY_TYPE_MAPPINGS_HPP
#define OBJECTMATRIX_PEER_NUMPY_TYPE_MAPPINGS_HPP

#include <numpy/arrayobject.h>

namespace peer {

namespace util
{
    namespace detail
    {
        namespace for_compile_time_error
        {
            template<typename T_>
            class numpy_does_not_support_the_type_;
        }
    }

    template<typename T_>
    struct get_numpy_typecode {
        static const std::size_t _ = sizeof(
                detail::for_compile_time_error::
                numpy_does_not_support_the_type_<T_>);
    };

#   define DEFINE_NUMPY_TYPECODE_ASSOC(__type__, __value__)  \
        template<> struct get_numpy_typecode<__type__> \
        { \
            BOOST_STATIC_CONSTANT(enum NPY_TYPES, value = __value__); \
        }
        DEFINE_NUMPY_TYPECODE_ASSOC(bool,            NPY_BOOL);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_byte,        NPY_BYTE);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_ubyte,       NPY_UBYTE);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_short,       NPY_SHORT);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_ushort,      NPY_USHORT);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_int,         NPY_INT);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_uint,        NPY_UINT);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_long,        NPY_LONG);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_ulong,       NPY_ULONG);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_longlong,    NPY_LONGLONG);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_ulonglong,   NPY_ULONGLONG);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_float,       NPY_FLOAT);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_double,      NPY_DOUBLE);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_longdouble,  NPY_LONGDOUBLE);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_cfloat,      NPY_CFLOAT);
        DEFINE_NUMPY_TYPECODE_ASSOC(std::complex<npy_float>, NPY_CFLOAT);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_cdouble,     NPY_CDOUBLE);
        DEFINE_NUMPY_TYPECODE_ASSOC(std::complex<npy_double>, NPY_CDOUBLE);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_clongdouble, NPY_CLONGDOUBLE);
        DEFINE_NUMPY_TYPECODE_ASSOC(
            std::complex<npy_longdouble>, NPY_CLONGDOUBLE);
        DEFINE_NUMPY_TYPECODE_ASSOC(boost::python::object, NPY_OBJECT);
        DEFINE_NUMPY_TYPECODE_ASSOC(std::string,           NPY_STRING);
#   define TMP std::basic_string<wchar_t, std::char_traits<wchar_t> >
        DEFINE_NUMPY_TYPECODE_ASSOC(TMP, NPY_UNICODE);
#   undef TMP
        DEFINE_NUMPY_TYPECODE_ASSOC(void,                  NPY_VOID);
        DEFINE_NUMPY_TYPECODE_ASSOC(char,                  NPY_CHAR);
#   undef DEFINE_NUMPY_TYPECODE_ASSOC
} // namespace util

} // namespace peer

#endif /* OBJECTMATRIX_PEER_NUMPY_TYPE_MAPPINGS_HPP */
