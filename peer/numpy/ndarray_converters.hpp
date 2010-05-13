#ifndef PEER_NUMPY_NDARRAY_CONVERTERS_HPP
#define PEER_NUMPY_NDARRAY_CONVERTERS_HPP

#include <stdexcept>
#include <complex>
#include <vector>
#include <boost/python.hpp>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <numpy/arrayobject.h>
#include "peer/utils.hpp"
#include "peer/numpy/pyarray_backed_allocator.hpp"
#include "peer/numpy/type_mappings.hpp"

namespace peer {

namespace util
{
    namespace detail
    {
        namespace for_compile_time_error
        {
            template<typename T_>
            class numpy_type_does_not_yield_from_;
        };

        template<typename Tarray_>
        struct to_ndarray_converter
        {
            typedef for_compile_time_error::
                    numpy_type_does_not_yield_from_<Tarray_> _;
        };

        template<typename T_, std::size_t N_>
        struct to_ndarray_converter<
                boost::multi_array<T_, N_, pyarray_backed_allocator<T_> > >
        {
            typedef boost::multi_array<T_, N_, pyarray_backed_allocator<T_> > source_type;
            static PyObject* convert(const source_type& val)
            {
                const npy_intp* dims;
                npy_intp _dims[N_];

                if (sizeof(npy_intp) == sizeof(typename source_type::size_type))
                {
                    dims = reinterpret_cast<const npy_intp*>(val.shape());
                }
                else
                {
                    for (std::size_t i = 0; i < N_; ++i)
                        _dims[i] = val.shape()[i];
                    dims = _dims;
                }
                PyObject* retval = PyArray_New(&PyArray_Type, N_,
                        const_cast<npy_intp*>(dims),
                        get_numpy_typecode<T_>::value, NULL,
                        const_cast<source_type&>(val).origin(), 0,
                        NPY_CARRAY, NULL);
                reinterpret_cast<PyArrayObject*>(retval)->flags |= NPY_OWNDATA;
                return retval;
            }
        };

        template<typename T_>
        struct to_ndarray_converter<std::vector<T_, pyarray_backed_allocator<T_ > > >
        {
            typedef std::vector<T_, pyarray_backed_allocator<T_> > source_type;

            static PyObject* convert(const source_type& val)
            {
                const npy_intp dims[1] = { val.size() };

                PyObject* retval = PyArray_New(&PyArray_Type, 1,
                        const_cast<npy_intp*>(dims),
                        get_numpy_typecode<T_>::value, NULL,
                        //const_cast<source_type&>(val).data(), 0,
                        &const_cast<source_type&>(val)[0], 0,
                        NPY_CARRAY, NULL);
                reinterpret_cast<PyArrayObject*>(retval)->flags |= NPY_OWNDATA;
                return retval;
            }
        };


        template<typename T_, typename Talloc_>
        struct to_ndarray_converter<std::vector<T_, Talloc_> >
        {
            typedef std::vector<T_, Talloc_> source_type;

            static PyObject* convert(const source_type& val)
            {
                typedef pyarray_backed_allocator<boost::python::object>
                        pyobject_array_allocator_type;

                BOOST_STATIC_ASSERT(
                        sizeof(boost::python::object) == sizeof(PyObject*)); 

                const npy_intp dims[1] = { val.size() };
                pyobject_array_allocator_type alloc(false);
                boost::python::object* converted_data =
                        alloc.allocate(val.size());

                boost::python::object* di = converted_data; 
                try
                {
                    for (typename source_type::const_iterator i(val.begin()),
                                                              e(val.end());
                            i != e; ++i)
                    {
                        new(di) boost::python::object(*i);
                        ++di; // this must be incremented after the pointed
                              // object pointed is successfully initialized
                    }
                }
                catch (const std::exception&)
                {
                    std::for_each(converted_data, di, functor::destructor_fun<
                            pyobject_array_allocator_type>(alloc));
                    return NULL; 
                }

                PyObject* retval = PyArray_New(&PyArray_Type, 1,
                        const_cast<npy_intp*>(dims),
                        get_numpy_typecode<boost::python::object>::value,
                        NULL, converted_data, 0,
                        NPY_CARRAY, NULL);
                 reinterpret_cast<PyArrayObject*>(retval)->flags |= NPY_OWNDATA;
                 return retval;
            }
        };

        template<typename T_, typename Telem_, std::size_t N_>
        struct array_to_ndarray_converter
        {
            typedef T_ native_type;
            
            static PyObject* convert( const native_type& p )
            {
                static const npy_intp dims[1] = { N_ };
                void* data(PyDataMem_NEW(N_ * sizeof(Telem_)));
                memcpy(data, static_cast<const void*>( &p[0] ),
                       N_ * sizeof(Telem_));
                PyObject* array( PyArray_New(&PyArray_Type, 1, 
                                             const_cast<npy_intp*>(dims),
                                             peer::util::get_numpy_typecode<
                                                 Telem_>::value,
                                             NULL, data, 0, NPY_CARRAY, NULL));
                reinterpret_cast<PyArrayObject*>(array)->flags |= NPY_OWNDATA;
                return array;
            }
        };

        template<typename T_, std::size_t N_>
        struct to_ndarray_converter<boost::array<T_, N_> >
            : array_to_ndarray_converter<boost::array<T_, N_>, T_, N_> {};

        template<typename T_, std::size_t N_>
        struct to_ndarray_converter<T_[N_]>
            : array_to_ndarray_converter<T_[N_], T_, N_> {};

    } // namespace detail

    template<typename Tarray_>
    inline void register_multi_array_converter()
    {
        static bool registered = false;
        if (!registered)
        {
            boost::python::to_python_converter<
                Tarray_, detail::to_ndarray_converter<Tarray_> >();
            registered = true;
        }
    }
} // namespace util

} // namespace peer

#endif /* PPER_NUMPY_NDARRAY_CONVERTERS_HPP */
