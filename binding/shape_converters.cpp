#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/python.hpp>
#include <boost/tuple/tuple.hpp>
#include "peer/utils.hpp"

#include "binding_common.hpp"

namespace binding {

struct sphere_to_python_converter
{
    typedef Sphere native_type;

    static PyObject* convert(native_type const& v)
    {
        return boost::python::incref(
            boost::python::object(boost::make_tuple(
                v.position(), v.radius())).ptr());
    }
};

struct python_to_sphere_converter
{
    typedef Sphere native_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PyTuple_Check(pyo))
        {
            return 0;
        }
        if (PyTuple_Size(pyo) != 2)
        {
            return 0;
        }
        return pyo;
    }

    static void construct(PyObject* pyo, 
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        PyObject* items[] = { PyTuple_GetItem(pyo, 0), PyTuple_GetItem(pyo, 1) };
        void* storage(reinterpret_cast<
            boost::python::converter::rvalue_from_python_storage<
                native_type >*
            >(data)->storage.bytes);
        new (storage) native_type (
            boost::python::extract<
                native_type::position_type>(items[0]),
            PyFloat_AsDouble(items[1]));
        data->convertible = storage;
    }
};

void register_sphere_converters()
{
    peer::util::to_native_converter<Sphere, python_to_sphere_converter>();
}

} // namespace binding
