#ifndef OBJECTMATRIX_PEER_PICKLE_SUPPORT
#define OBJECTMATRIX_PEER_PICKLE_SUPPORT

#include <Python.h>
#include <pyerrors.h>
#include <import.h>
#include <boost/python.hpp>

namespace peer { namespace pickle {

static const char reconstruct_func_name[] = "__reconstruct__";

static PyObject* reconstruct(PyObject* self, PyObject* args)
{
    PyTypeObject* klass;
    PyObject* base;
    PyObject* state;

    if (!PyArg_ParseTuple(args, "OOO", &klass, &base, &state))
        return NULL;

    if (!PyType_Check(klass)) {
        PyErr_SetString(PyExc_TypeError, "argument 1 must be a type object");
        return NULL;
    }

    if (!PyType_Check(base)) {
        PyErr_SetString(PyExc_TypeError, "argument 2 must be a type object");
        return NULL;
    }

    if (!PyTuple_Check(state)) {
        PyErr_SetString(PyExc_TypeError, "argument 3 must be a tuple");
        return NULL;
    }

    return klass->tp_new(klass, state, NULL);
}


static PyObject* reduce(PyObject* self) try
{
    using namespace boost::python;

    BOOST_ASSERT(self->ob_type);
    BOOST_ASSERT(self->ob_type->tp_base);

    BOOST_ASSERT(self != Py_None);

    object state(getattr(object(borrowed(self)), "__getstate__")());
    object module(borrowed(PyImport_Import(
                getattr(object(borrowed(
                    reinterpret_cast<PyObject*>(self->ob_type))),
                    "__module__").ptr())));

    return incref(make_tuple(
        getattr(module, reconstruct_func_name),
        make_tuple(
            borrowed(incref(reinterpret_cast<PyObject*>(self->ob_type))),
            borrowed(incref(reinterpret_cast<PyObject*>(self->ob_type->tp_base))),
            state)).ptr());
}
catch (boost::python::error_already_set const&)
{
    return NULL;
}

static inline void register_reconstructor()
{
    using namespace boost::python;
    static bool registered = false;
    if (!registered) {
        static PyMethodDef def = {
            const_cast<char*>(reconstruct_func_name),
            &reconstruct,
            METH_VARARGS, const_cast<char*>("")
        };
        scope().attr(reconstruct_func_name) = borrowed(PyCFunction_NewEx(
            &def, NULL, getattr(scope(), "__name__").ptr()));
    }
}

} } // namespace peer::pickle

#endif /* OBJECTMATRIX_PEER_PICKLE_SUPPORT */
