#ifndef OBJECTMATRIX_PEER_SPHERE_HPP
#define OBJECTMATRIX_PEER_SPHERE_HPP

#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#include "sphere.hpp"
#include "peer/utils.hpp"
#include "peer/numpy/type_mappings.hpp"

namespace peer {

class Sphere
{
public:
    typedef ::sphere<double> impl_type;
    typedef impl_type::value_type value_type;

public:
    Sphere(): impl_(impl_type::position_type(0, 0, 0), 0) {}

    Sphere(const impl_type::position_type& p, value_type r): impl_(p, r) {}

    PyObject* __repr__() const
    {
        return util::pystr_from_repr(&impl_);
    }

    value_type _get_x() const
    {
        return impl_.position.x();
    }

    void _set_x(value_type val)
    {
        impl_.position.x() = val;
    }

    value_type _get_y() const
    {
        return impl_.position.y();
    }

    void _set_y(value_type val)
    {
        impl_.position.y() = val;
    }

    value_type _get_z() const
    {
        return impl_.position.z();
    }

    void _set_z(value_type val)
    {
        impl_.position.z() = val;
    }

    value_type _get_radius() const
    {
        return impl_.radius;
    }

    void _set_radius(value_type val)
    {
        impl_.radius = val;
    }

    operator impl_type&()
    {
        return impl_;
    }

    operator const impl_type&() const
    {
        return impl_;
    }

    struct position_to_ndarray_converter
    {
        typedef impl_type::position_type native_type;

        static PyObject* convert(const native_type& p)
        {
            static const npy_intp dims[1] = { native_type::size() };
            return PyArray_New(&PyArray_Type, 1, const_cast<npy_intp*>(dims),
                   util::get_numpy_typecode<value_type>::value, NULL,
                    const_cast<void*>(static_cast<const void*>(&p)),
                    0, NPY_CARRAY, NULL);
        }
    };

    struct ndarray_to_position_converter
    {
        typedef impl_type::position_type native_type;

        static void* convertible(PyObject* ptr)
        {
            if (!PyArray_Check(ptr))
            {
                return NULL;
            }

            PyObject* retval(PyArray_CastToType(
                    reinterpret_cast<PyArrayObject*>(ptr),
                    PyArray_DescrFromType(
                        util::get_numpy_typecode<
                            native_type::value_type>::value), 0));
            if (!retval)
            {
                return NULL;
            }

            if (PyObject_Size(retval) != 3)
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
            data->stage1.convertible = new(data->storage.bytes) native_type(
                    reinterpret_cast<double*>(PyArray_DATA(array_obj)));
            boost::python::decref(reinterpret_cast<PyObject*>(array_obj));
        }
    };

    struct seq_to_position_converter
    {
        typedef impl_type::position_type native_type;

        static void* convertible(PyObject* ptr)
        {
            if (!PySequence_Check(ptr))
            {
                return NULL;
            }

            if (PySequence_Size(ptr) != 3)
            {
                return NULL;
            }

            return ptr;
        }

        static void construct(PyObject* ptr,
                boost::python::converter::rvalue_from_python_storage<native_type>* data)
        {
            data->stage1.convertible = new(data->storage.bytes) native_type(
                PyFloat_AsDouble(PySequence_GetItem(ptr, 0)),
                PyFloat_AsDouble(PySequence_GetItem(ptr, 1)),
                PyFloat_AsDouble(PySequence_GetItem(ptr, 2)));
        }
    };

    inline static void __register_class()
    {
        using namespace boost::python;

        to_python_converter<impl_type::position_type,
                position_to_ndarray_converter>();
        util::to_native_converter<impl_type::position_type,
                ndarray_to_position_converter>();
        util::to_native_converter<impl_type::position_type,
                seq_to_position_converter>();

        class_<Sphere>("Sphere")
            .def(init<const ::position<double>&, double>())
            .def("__repr__", &Sphere::__repr__)
            .add_property("x", &Sphere::_get_x, &Sphere::_set_x)
            .add_property("y", &Sphere::_get_y, &Sphere::_set_y)
            .add_property("z", &Sphere::_get_z, &Sphere::_set_z)
            .add_property("radius", &Sphere::_get_radius, &Sphere::_set_radius);
    }
private:
    impl_type impl_;
};

} // namespace peer
#endif /* OBJECTMATRIX_PEER_SPHERE_HPP */
