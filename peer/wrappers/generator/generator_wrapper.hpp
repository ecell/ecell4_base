#ifndef PEER_WRAPPERS_GENERATOR_GENERATOR_WRAPPER_HPP
#define PEER_WRAPPERS_GENERATOR_GENERATOR_WRAPPER_HPP

#include <Python.h>
#include "generator.hpp"
#include "peer/utils.hpp"

namespace peer { namespace wrappers {

template<typename Tgen_, typename Trcg_ = boost::python::return_by_value >
class generator_wrapper
{
public:
    typedef Tgen_ generator_type;
    typedef typename Trcg_::template apply<typename generator_type::result_type>::type result_converter_type;

public:
    void* operator new(size_t)
    {
        PyObject* retval = PyObject_New(PyObject, &__class__);
        return retval;
    }

    void operator delete(void* ptr)
    {
        reinterpret_cast<PyObject*>(ptr)->ob_type->tp_free(reinterpret_cast<PyObject*>(ptr));
    }

    Tgen_ const& ptr() const
    {
        return impl_;
    }

    Tgen_& ptr()
    {
        return impl_;
    }

    generator_wrapper(Tgen_ const& impl): impl_(impl) {}

    ~generator_wrapper() {}

    static PyObject* create(Tgen_ const& impl)
    {
        return reinterpret_cast<PyObject*>(new generator_wrapper(impl));
    } 

    static void __dealloc__(generator_wrapper* self)
    {
        delete self;
    }

    static PyObject* next(generator_wrapper* self)
    {
        using namespace boost::python;

        if (!valid(self->impl_))
        {
            PyErr_SetNone(PyExc_StopIteration);
            return NULL;
        }
        return result_converter_type()(self->impl_());
    }


    static PyObject* __length_hint__(generator_wrapper* self)
    {
        try
        {
#ifdef HAVE_PYINT_FROMSIZE_T
            return PyInt_FromSize_t(count(self->impl_));
#else
            std::size_t i(count(self->impl_));
            if (i >= static_cast<std::size_t>(LONG_MIN) &&
                    i <= static_cast<std::size_t>(LONG_MAX))
            {
                return PyInt_FromLong((long)i);
            }
            return PyLong_FromUnsignedLongLong(i);
#endif
        }
        catch (std::exception const&) {}

        PyErr_SetNone(PyExc_TypeError);
        return NULL;
    }

    static PyTypeObject* __class_init__(const char* name, PyObject* mod)
    {
        using namespace boost::python;
        if (__name__.empty())
        {
            __name__ = mod && PyModule_Check(mod) ?
                extract<std::string>(object(borrowed(mod)).attr("__name__"))()
                + "." + name: std::string(name);
            __class__.tp_name = const_cast<char*>(__name__.c_str());
            PyType_Ready(&__class__);
        }
        return &__class__;
    }

    static void __register_class(const char* name);

protected:
    PyObject_VAR_HEAD
    static PyTypeObject __class__;
    static PyMethodDef __methods__[];
    static std::string __name__;
    Tgen_ impl_;
};

template<typename Tgen_, typename Trcg_>
inline void generator_wrapper<Tgen_, Trcg_>::__register_class(const char* name)
{
    using namespace boost::python;
    PyTypeObject* klass(generator_wrapper::__class_init__(name, reinterpret_cast<PyObject*>(scope().ptr())));
    Py_INCREF(klass);
    scope().attr(name) = object(borrowed(reinterpret_cast<PyObject*>(klass)));
}

template<typename Tgen_, typename Trcg_>
std::string generator_wrapper<Tgen_, Trcg_>::__name__;

template<typename Tgen_, typename Trcg_>
PyMethodDef generator_wrapper<Tgen_, Trcg_>::__methods__[] = {
    { "__length_hint__", (PyCFunction)generator_wrapper::__length_hint__, METH_NOARGS, "" },
    { NULL, NULL }
};

template<typename Tgen_, typename Trcg_>
PyTypeObject generator_wrapper<Tgen_, Trcg_>::__class__ = {
    PyObject_HEAD_INIT(&PyType_Type)
    0,                  /* ob_size */
    0,                  /* tp_name */
    sizeof(generator_wrapper), /* tp_basicsize */
    0,                  /* tp_itemsize */
    /* methods */
    (destructor)&generator_wrapper::__dealloc__, /* tp_dealloc */
    0,                  /* tp_print */
    0,                  /* tp_getattr */
    0,                  /* tp_setattr */
    0,                  /* tp_compare */
    0,                  /* tp_repr */
    0,                  /* tp_as_number */
    0,                  /* tp_as_sequence */
    0,                  /* tp_as_mapping */
    0,                  /* tp_hash */
    0,                  /* tp_call */
    0,                  /* tp_str */
    PyObject_GenericGetAttr,        /* tp_getattro */
    0,                  /* tp_setattro */
    0,                  /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER, /* tp_flags */
    0,                  /* tp_doc */
    0,                  /* tp_traverse */
    0,                  /* tp_clear */
    0,                  /* tp_richcompare */
    0,                  /* tp_weaklistoffset */
    PyObject_SelfIter,  /* tp_iter */
    (iternextfunc)&generator_wrapper::next,  /* tp_iternext */
    generator_wrapper::__methods__,  /* tp_methods */
    0,                  /* tp_members */
    0,                  /* tp_getset */
    0,                  /* tp_base */
    0,                  /* tp_dict */
    0,                  /* tp_descr_get */
    0,                  /* tp_descr_set */
    0,                  /* tp_dictoffset */
    0,                  /* tp_init */
    0,                  /* tp_alloc */
    0,                  /* tp_new */
    0                   /* tp_free */
};

} } // namespace peer::util

#endif /* PEER_WRAPPERS_GENERATOR_GENERATOR_WRAPPER_HPP */
