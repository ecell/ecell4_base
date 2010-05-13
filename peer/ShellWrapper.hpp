#ifndef PEER_SHELL_WRAPPER_HPP
#define PEER_SHELL_WRAPPER_HPP

#include <cstddef>
#include <string>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include <boost/python.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/object/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <numpy/arrayobject.h>

#include <unistd.h>

#include "numpy/type_mappings.hpp"
#include "pickle_support.hpp"

namespace peer {

template<typename Timpl_>
class ShellWrapper
{
public:
    struct to_python_converter
    {
        static PyObject* convert(Timpl_ const& impl)
        {
            return ShellWrapper::create(impl);
        }

        static PyTypeObject* get_pytype()
        {
            return &ShellWrapper::__class__;
        }
    };

    struct to_native_converter
    {
        static void* convertible(PyObject* pyo)
        {
            if (!PyObject_TypeCheck(pyo, &ShellWrapper::__class__))
            {
                return 0;
            }
            return pyo;
        }

        static void construct(PyObject* pyo,
                              boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            void* storage(reinterpret_cast<
                boost::python::converter::rvalue_from_python_storage<Timpl_>* >(
                    data)->storage.bytes);
            new (storage) Timpl_(reinterpret_cast<ShellWrapper*>(pyo)->impl_);
            data->convertible = storage;
        }
    };

public:
    static PyTypeObject* __class_init__(const char* name, PyObject* mod)
    {
        using namespace boost::python;
        __name__ = static_cast<std::string>(
            extract<std::string>(object(borrowed(mod)).attr("__name__")))
            + "." + name;
        __class__.tp_name = const_cast<char*>(__name__.c_str());
        PyType_Ready(&__class__);
        dict(borrowed(__class__.tp_dict))["__safe_for_unpickling__"] = object(true);
        return &__class__;
    }

    static PyObject* create()
    {
        return reinterpret_cast<PyObject*>(new ShellWrapper());
    }

    static PyObject* create(Timpl_ const& impl)
    {
        return reinterpret_cast<PyObject*>(new ShellWrapper(impl));
    }

    static PyObject* __new__(PyTypeObject* klass, PyObject* arg, PyObject* kwarg)
    {
        PyObject* retval = NULL;
        switch (PyTuple_Size(arg))
        {
        default:
            PyErr_SetString(PyExc_TypeError, "the number of arguments must be either 0 or 2");
            return NULL;

        case 2:
            retval = create();
            if (set_did(reinterpret_cast<ShellWrapper*>(retval), 
                        PyTuple_GetItem(arg, 0), 0)
                || set_shape(reinterpret_cast<ShellWrapper*>(retval), 
                        PyTuple_GetItem(arg, 1), 0))
            {
                ShellWrapper::operator delete(retval);
                return NULL;
            }
            break;

        case 0:
            retval = create();
            break;
        }

        if (PyErr_Occurred())
        {
            boost::python::decref(reinterpret_cast<PyObject*>(retval));
            return NULL;
        }
        return reinterpret_cast<PyObject*>(retval);
    }

    static PyObject* __str__(ShellWrapper* self)
    {
        std::string retval(boost::lexical_cast<std::string>(self->impl_));
        return PyString_FromStringAndSize(retval.data(), retval.size());
    }

    static PyObject* __repr__(ShellWrapper* self)
    {
        return __str__(self);
    }

    static void __dealloc__(ShellWrapper* self)
    {
        delete self;
    }

    static PyObject* get_did(ShellWrapper* self)
    {
        return boost::python::incref(
            boost::python::object(self->impl_.did()).ptr());
    }

    static int set_did(ShellWrapper* self, PyObject* val, void *)
    try
    {
        self->impl_.did() = boost::python::extract<typename Timpl_::domain_id_type>(val);
        return 0;
    }
    catch (boost::python::error_already_set const&)
    {
        return -1;
    }

    static PyObject* get_shape(ShellWrapper* self)
    {
        return boost::python::incref(boost::python::object(self->impl_.shape()).ptr());
    }

    static int set_shape(ShellWrapper* self, PyObject* val, void *)
    try
    {
        self->impl_.shape() = boost::python::extract<typename Timpl_::shape_type>(val)();
        return 0;
    }
    catch (boost::python::error_already_set const&)
    {
        return -1;
    }

    static PyObject* __getstate__(ShellWrapper* self)
    try
    {
        return boost::python::incref(
            boost::python::make_tuple(
                boost::python::borrowed(get_did(self)),
                boost::python::borrowed(get_shape(self))).ptr());
    }
    catch (boost::python::error_already_set const&)
    {
        return NULL;
    }

    static PyObject* __reduce__(ShellWrapper* self)
    {
        return pickle::reduce(reinterpret_cast<PyObject*>(self));
    }

    static PyObject* __reduce_ex__(ShellWrapper* self, PyObject* arg)
    {
        return pickle::reduce(reinterpret_cast<PyObject*>(self));
    }

    static long __hash__(ShellWrapper* self)
    {
#if defined(HAVE_TR1_FUNCTIONAL)
        using namespace std::tr1;
#elif defined(HAVE_STD_HASH)
        using namespace std;
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
        using namespace boost;
#endif
        return static_cast<long>(hash<Timpl_>()(self->impl_));
    }

    // this method borrows a reference to "name" from the caller.
    // i.e. name should be statically allocated
    static void __register_class(char const* name)
    {
        using namespace boost::python;
        pickle::register_reconstructor();
        PyTypeObject* klass(ShellWrapper::__class_init__(name, reinterpret_cast<PyObject*>(scope().ptr())));
        Py_INCREF(klass);
        scope().attr(name) = object(borrowed(reinterpret_cast<PyObject*>(klass)));
        util::to_native_converter<Timpl_, to_native_converter>();
        boost::python::to_python_converter<Timpl_, to_python_converter>();
    }

    static PyTypeObject* get_class()
    {
        return &__class__;
    }

    Timpl_& wrapped()
    {
        return impl_;
    }

    Timpl_ const& wrapped() const
    {
        return impl_;
    }

protected:
    ShellWrapper(): impl_() {}

    ShellWrapper(Timpl_ const& impl): impl_(impl) {}

    void* operator new(size_t)
    {
        PyObject* retval = PyObject_New(PyObject, &__class__);
        return retval;
    }

    void operator delete(void* ptr)
    {
        reinterpret_cast<PyObject*>(ptr)->ob_type->tp_free(reinterpret_cast<PyObject*>(ptr));
    }

    ~ShellWrapper() {}

protected:
    PyObject_VAR_HEAD
    static PyTypeObject __class__;
    static PyGetSetDef __getsets__[];
    static PyMethodDef __methods__[];
    static PySequenceMethods __sequence_methods__;
    static boost::python::object reconstructor;
    static std::string __name__;
    Timpl_ impl_;
};

template<typename Timpl_>
std::string ShellWrapper<Timpl_>::__name__;

template<typename Timpl_>
PyMethodDef ShellWrapper<Timpl_>::__methods__[] = {
    { "__getstate__", (PyCFunction)ShellWrapper::__getstate__, METH_NOARGS, "" },
    { "__reduce__", (PyCFunction)ShellWrapper::__reduce__, METH_NOARGS, "" },
    { "__reduce_ex__", (PyCFunction)ShellWrapper::__reduce_ex__, METH_O, "" },
    { NULL, NULL }
};

template<typename Timpl_>
PyGetSetDef ShellWrapper<Timpl_>::__getsets__[] = {
    {
        const_cast<char*>("shape"),
        (getter)ShellWrapper::get_shape,
        (setter)ShellWrapper::set_shape,
        const_cast<char*>("")
    },
    {
        const_cast<char*>("did"),
        (getter)ShellWrapper::get_did,
        (setter)ShellWrapper::set_did,
        const_cast<char*>("")
    },
    { NULL }
};

template<typename Timpl_>
PyTypeObject ShellWrapper<Timpl_>::__class__ = {
    PyObject_HEAD_INIT(&PyType_Type)
    0,                  /* ob_size */
    0,                  /* tp_name */
    sizeof(ShellWrapper), /* tp_basicsize */
    0,                  /* tp_itemsize */
    /* methods */
    (destructor)&ShellWrapper::__dealloc__, /* tp_dealloc */
    0,                  /* tp_print */
    0,                  /* tp_getattr */
    0,                  /* tp_setattr */
    0,                  /* tp_compare */
    (reprfunc)&ShellWrapper::__repr__,                    /* tp_repr */
    0,                  /* tp_as_number */
    0,                  /* tp_as_sequence */
    0,                  /* tp_as_mapping */
    (hashfunc)&ShellWrapper::__hash__,                    /* tp_hash */
    0,                  /* tp_call */
    (reprfunc)&ShellWrapper::__str__,                    /* tp_str */
    PyObject_GenericGetAttr,        /* tp_getattro */
    0,                  /* tp_setattro */
    0,                  /* tp_as_buffer */
    Py_TPFLAGS_HAVE_CLASS | Py_TPFLAGS_HAVE_RICHCOMPARE,/* tp_flags */
    0,                  /* tp_doc */
    0,                  /* tp_traverse */
    0,                  /* tp_clear */
    0,                  /* tp_richcompare */
    0,                  /* tp_weaklistoffset */
    0,                  /* tp_iter */
    0,                  /* tp_iternext */
    ShellWrapper::__methods__,                    /* tp_methods */
    0,                  /* tp_members */
    ShellWrapper::__getsets__, /* tp_getset */
    &PyBaseObject_Type, /* tp_base */
    0,                  /* tp_dict */
    0,                  /* tp_descr_get */
    0,                  /* tp_descr_set */
    0,                  /* tp_dictoffset */
    0,                  /* tp_init */
    0,                  /* tp_alloc */
    ShellWrapper::__new__,  /*tp_new */
    0                   /* tp_free */
};

} //namespace peer

#endif /* PEER_SHELL_WRAPPER_HPP */
