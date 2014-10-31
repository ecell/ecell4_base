#ifndef WRAPPED_MULTI_ARRAY_HPP
#define WRAPPED_MULTI_ARRAY_HPP

#include <boost/multi_array.hpp>
#include <boost/intrusive_ptr.hpp>

#include "peer/numpy/type_mappings.hpp"

namespace peer { namespace util {

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


template<typename T_, std::size_t Ndims_>
class ndarray_wrapped_multi_array_converter
{
public:
    typedef python_array_lifecycle_manager<T_> lcmgr_type;
    typedef wrapped_multi_array<T_, Ndims_> native_type;

public:
    static void* convertible(PyObject* ptr)
    {
        if (!PyArray_Check(ptr))
        {
            return NULL;
        }

        PyObject* retval(
            PyArray_CastToType(
                reinterpret_cast<PyArrayObject*>(ptr),
                PyArray_DescrFromType(
                    get_numpy_typecode<
                        typename native_type::element>::value), 0));
        if (!retval)
        {
            return NULL;
        }

        if (PyArray_NDIM(reinterpret_cast<PyArrayObject*>(retval)) != Ndims_)
        {
            boost::python::decref(retval);
            return NULL;
        }

        return retval;
    }

    static void construct(PyObject* ptr,
            boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        PyArrayObject* array_obj = static_cast<PyArrayObject*>(
                data->stage1.convertible);
        typename native_type::index_list ma_strides;

        for (std::size_t i = 0; i < Ndims_; ++i)
        {
            ma_strides[i] = array_obj->strides[i] / sizeof(T_);
        }

        data->stage1.convertible = new(data->storage.bytes) native_type(
                new lcmgr_type(array_obj),
                *reinterpret_cast<boost::array<npy_intp, Ndims_>*>(
                    array_obj->dimensions),
                ma_strides,
                PyArray_ISCONTIGUOUS(array_obj) ?
                    static_cast<boost::general_storage_order<Ndims_> >(
                        boost::c_storage_order()):
                    static_cast<boost::general_storage_order<Ndims_> >(
                        boost::fortran_storage_order())
                );
    }
};


template<typename T_>
class seq_wrapped_multi_array_converter
{
public:
    typedef python_array_lifecycle_manager<T_> lcmgr_type;
    typedef wrapped_multi_array<T_, 1> native_type;

public:
    static void* convertible(PyObject* ptr)
    {
        if (!PySequence_Check(ptr))
        {
            return NULL;
        }

        return ptr;
    }

    static void construct(PyObject* ptr,
            boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        PyObject* array_obj = PyArray_FromObject(
            static_cast<PyObject*>(data->stage1.convertible),
            get_numpy_typecode<T_>::value, 1, 1);
        static typename native_type::index_list ma_strides = { { 1 } };

        data->stage1.convertible = new(data->storage.bytes) native_type(
                new lcmgr_type(reinterpret_cast<PyArrayObject *>(array_obj)),
                *reinterpret_cast<boost::array<npy_intp, 1>*>(
                    PyArray_DIMS(array_obj)),
                ma_strides,
                static_cast<boost::general_storage_order<1> >(
                        boost::c_storage_order()));
    }
};


template<typename T_, std::size_t Ndims_>
void register_ndarray_wrapped_multi_array_converter()
{
    typedef ndarray_wrapped_multi_array_converter<T_, Ndims_> Converter;
    boost::python::converter::registry::push_back(
        &Converter::convertible,
        reinterpret_cast<boost::python::converter::constructor_function>(
            &Converter::construct),
        boost::python::type_id<typename Converter::native_type>());
}


template<typename T_>
void register_seq_wrapped_multi_array_converter()
{
    typedef seq_wrapped_multi_array_converter<T_> Converter;
    boost::python::converter::registry::push_back(
        &Converter::convertible,
        reinterpret_cast<boost::python::converter::constructor_function>(
            &Converter::construct),
        boost::python::type_id<typename Converter::native_type>());
}

} } // namespace peer::util

#endif
