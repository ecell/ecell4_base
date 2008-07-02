#ifndef WRAPPED_MULTI_ARRAY_HPP
#define WRAPPED_MULTI_ARRAY_HPP

#include <boost/multi_array.hpp>
#include <boost/intrusive_ptr.hpp>


namespace for_compile_time_error
{
template<typename T_>
class numpy_does_not_support_the_type_;
}





template<typename T_>
struct get_numpy_typecode {
    static const std::size_t value = sizeof(
            for_compile_time_error::
            numpy_does_not_support_the_type_<T_>);
};

#define DEFINE_NUMPY_TYPECODE_ASSOC(__type__, __value__)  \
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
#define TMP std::basic_string<wchar_t, std::char_traits<wchar_t> >
DEFINE_NUMPY_TYPECODE_ASSOC(TMP, NPY_UNICODE);
#undef TMP
DEFINE_NUMPY_TYPECODE_ASSOC(void,                  NPY_VOID);
DEFINE_NUMPY_TYPECODE_ASSOC(char,                  NPY_CHAR);
#undef DEFINE_NUMPY_TYPECODE_ASSOC


template<typename T_>
class lifecycle_manager
{
public:
    virtual ~lifecycle_manager() {}

    virtual T_* data() = 0;

    virtual const T_* data() const = 0;

    inline T_* operator*()
    {
        return data();
    }

    inline const T_* operator*() const
    {
        return data();
    }
};

template<typename T_>
class python_array_lifecycle_manager: public lifecycle_manager<T_>
{
private:
    typedef python_array_lifecycle_manager self_type;

public:
    python_array_lifecycle_manager(PyArrayObject* array_obj)
        : array_obj_(array_obj) {}

    python_array_lifecycle_manager(const self_type& that)
        : array_obj_(that.array_obj_)
    {
        Py_INCREF(array_obj_);
    }

    virtual ~python_array_lifecycle_manager()
    {
        Py_DECREF(array_obj_);
    }

    virtual T_* data()
    {
        return reinterpret_cast<T_*>(PyArray_DATA(array_obj_));
    }

    virtual const T_* data() const
    {
        return reinterpret_cast<T_*>(PyArray_DATA(array_obj_));
    }

private:
    PyArrayObject* array_obj_;
};

template<typename T, std::size_t NumDims>
class wrapped_multi_array: public boost::multi_array_ref<T, NumDims>
{
public:
    typedef boost::multi_array_ref<T, NumDims> super_type;
    typedef typename super_type::value_type value_type;
    typedef typename super_type::reference reference;
    typedef typename super_type::const_reference const_reference;
    typedef typename super_type::iterator iterator;
    typedef typename super_type::const_iterator const_iterator;
    typedef typename super_type::reverse_iterator reverse_iterator;
    typedef typename super_type::const_reverse_iterator const_reverse_iterator;
    typedef typename super_type::element element;
    typedef typename super_type::size_type size_type;
    typedef typename super_type::difference_type difference_type;
    typedef typename super_type::index index;
    typedef typename super_type::extent_range extent_range;
    typedef typename super_type::index_list index_list;
    typedef typename super_type::size_list size_list;
    typedef lifecycle_manager<T> lifecycle_manager_type;

    template <std::size_t NDims>
    struct const_array_view {
        typedef boost::detail::multi_array::const_multi_array_view<T,NDims> type;
    };

    template <std::size_t NDims>
    struct array_view {
        typedef boost::detail::multi_array::multi_array_view<T,NDims> type;
    };

    template <class ExtentList, class StrideList>
    wrapped_multi_array(
        lifecycle_manager_type* lm,
        ExtentList const& extents,
        StrideList const& strides,
        const boost::general_storage_order<NumDims>& so
        ): super_type(lm->data(), extents, so), lm_(lm)
    {
        boost::function_requires<
            boost::detail::multi_array::CollectionConcept<ExtentList> >();
        boost::function_requires<
            boost::detail::multi_array::CollectionConcept<StrideList> >();
        std::copy(strides.begin(), strides.end(),
            super_type::stride_list_.begin());
    }

private:
    boost::shared_ptr<lifecycle_manager_type> lm_;
};






#endif
