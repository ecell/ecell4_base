#ifndef BINDING_IDENTIFIER_HPP
#define BINDING_IDENTIFIER_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

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

#include <unistd.h>
#include <iostream>

#include "peer/pickle_support.hpp"
#include "peer/utils.hpp"

namespace binding {

template<typename Timpl_>
class IdentifierWrapper
{
public:
    struct to_python_converter
    {
        static PyObject* convert(Timpl_ const& impl)
        {
            return IdentifierWrapper::create(impl);
        }

        static PyTypeObject* get_pytype()
        {
            return &IdentifierWrapper::__class__;
        }
    };

    struct to_native_converter
    {
        static void* convertible(PyObject* pyo)
        {
            if (!PyObject_TypeCheck(pyo, &IdentifierWrapper::__class__))
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
            new (storage) Timpl_(reinterpret_cast<IdentifierWrapper*>(pyo)->impl_);
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

    static PyObject* create(Timpl_ const& impl)
    {
        return reinterpret_cast<PyObject*>(new IdentifierWrapper(impl));
    }

    static PyObject* __new__(PyTypeObject* klass, PyObject* arg, PyObject* kwarg)
    try
    {
        if (!PyTuple_Check(arg) || PyTuple_Size(arg) != 2)
        {
            PyErr_SetString(PyExc_TypeError, "the number of arguments must be 2");
            return NULL;
        }

        return create(Timpl_(
            typename Timpl_::value_type(
                boost::python::extract<typename Timpl_::lot_type>(PyTuple_GetItem(arg, 0)),
                boost::python::extract<typename Timpl_::serial_type>(PyTuple_GetItem(arg, 1))
            )));
    }
    catch (boost::python::error_already_set const&)
    {
        return NULL;
    }

    static PyObject* get_lot(IdentifierWrapper* self)
    {
        return boost::python::incref(boost::python::object(lot(self->impl_)).ptr());
    }

    static PyObject* get_serial(IdentifierWrapper* self)
    {
        return boost::python::incref(boost::python::object(serial(self->impl_)).ptr());
    }

    static PyObject* __str__(IdentifierWrapper* self)
    {
        std::string retval(boost::lexical_cast<std::string>(self->impl_));
        return PyString_FromStringAndSize(retval.data(), retval.size());
    }

    static PyObject* __repr__(IdentifierWrapper* self)
    {
        return __str__(self);
    }

    static void __dealloc__(IdentifierWrapper* self)
    {
        delete self;
    }

    static PyObject* __getstate__(IdentifierWrapper* self)
    try
    {
        return boost::python::incref(
            boost::python::make_tuple(
                boost::python::borrowed(get_lot(self)),
                boost::python::borrowed(get_serial(self))).ptr());
    }
    catch (boost::python::error_already_set const&)
    {
        return NULL;
    }

    static PyObject* __reduce__(IdentifierWrapper* self)
    {
        return peer::pickle::reduce(reinterpret_cast<PyObject*>(self));
    }

    static PyObject* __reduce_ex__(IdentifierWrapper* self, PyObject* arg)
    {
        return peer::pickle::reduce(reinterpret_cast<PyObject*>(self));
    }

    static long __hash__(IdentifierWrapper* self)
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

    static PyObject* __richcmp__(IdentifierWrapper* self, PyObject* rhs, int op)
    {
        PyObject* retval = Py_None;
        switch (op)
        {
        case Py_EQ:
            if (!PyObject_TypeCheck(rhs, &__class__))
            {
                retval = Py_False;
                Py_INCREF(retval);
            }
            else
            {
                IdentifierWrapper* const _rhs = reinterpret_cast<IdentifierWrapper*>(rhs);

                retval = PyBool_FromLong(self->impl_ == _rhs->impl_);
            }
            break;
        case Py_NE:
            if (!PyObject_TypeCheck(rhs, &__class__))
            {
                retval = Py_True;
                Py_INCREF(retval);
            }
            else
            {
                IdentifierWrapper* const _rhs = reinterpret_cast<IdentifierWrapper*>(rhs);
                retval = PyBool_FromLong(self->impl_ != _rhs->impl_);
            }
            break;
        default:
            if (!PyObject_TypeCheck(rhs, &__class__))
            {
                PyErr_SetString(PyExc_TypeError, "Comparison is not feasible");
                Py_RETURN_NONE;
            }
            IdentifierWrapper* const _rhs = reinterpret_cast<IdentifierWrapper*>(rhs);
            switch (op)
            {
            default:
                break;
            case Py_LT:
                retval = PyBool_FromLong(self->impl_ < _rhs->impl_);
                break;
            case Py_LE:
                retval = PyBool_FromLong(self->impl_ <= _rhs->impl_);
                break;
            case Py_GT:
                retval = PyBool_FromLong(self->impl_ > _rhs->impl_);
                break;
            case Py_GE:
                retval = PyBool_FromLong(self->impl_ >= _rhs->impl_);
                break;
            }
        }
        Py_INCREF(retval);
        return retval;
    }

    // this method borrows a reference to "name" from the caller.
    // i.e. name should be statically allocated
    static void __register_class(char const* name)
    {
        using namespace boost::python;
        peer::pickle::register_reconstructor();
        PyTypeObject* klass(IdentifierWrapper::__class_init__(name, reinterpret_cast<PyObject*>(scope().ptr())));
        Py_INCREF(klass);
        scope().attr(name) = object(borrowed(reinterpret_cast<PyObject*>(klass)));
        peer::util::to_native_converter<Timpl_, to_native_converter>();
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
    IdentifierWrapper(Timpl_ const& impl): impl_(impl) {}

    void* operator new(size_t)
    {
        PyObject* retval = PyObject_New(PyObject, &__class__);
        return retval;
    }

    void operator delete(void* ptr)
    {
        reinterpret_cast<PyObject*>(ptr)->ob_type->tp_free(reinterpret_cast<PyObject*>(ptr));
    }

    ~IdentifierWrapper() {}

protected:
    PyObject_VAR_HEAD
    static PyTypeObject __class__;
    static PyGetSetDef __getsets__[];
    static PyMethodDef __methods__[];
    static boost::python::object reconstructor;
    static std::string __name__;
    Timpl_ impl_;
};

template<typename Timpl_>
std::string IdentifierWrapper<Timpl_>::__name__;

template<typename Timpl_>
PyMethodDef IdentifierWrapper<Timpl_>::__methods__[] = {
    { "__getstate__", (PyCFunction)IdentifierWrapper::__getstate__, METH_NOARGS, "" },
    { "__reduce__", (PyCFunction)IdentifierWrapper::__reduce__, METH_NOARGS, "" },
    { "__reduce_ex__", (PyCFunction)IdentifierWrapper::__reduce_ex__, METH_O, "" },
    { NULL, NULL }
};

template<typename Timpl_>
PyGetSetDef IdentifierWrapper<Timpl_>::__getsets__[] = {
    {
        const_cast<char*>("lot"),
        (getter)IdentifierWrapper::get_lot,
        NULL,
        const_cast<char*>("")
    },
    {
        const_cast<char*>("serial"),
        (getter)IdentifierWrapper::get_serial,
        NULL,
        const_cast<char*>("")
    },
    { NULL }
};

template<typename Timpl_>
PyTypeObject IdentifierWrapper<Timpl_>::__class__ = {
	PyObject_HEAD_INIT(&PyType_Type)
	0,					/* ob_size */
	0,                  /* tp_name */
	sizeof(IdentifierWrapper), /* tp_basicsize */
	0,					/* tp_itemsize */
	/* methods */
	(destructor)&IdentifierWrapper::__dealloc__, /* tp_dealloc */
	0,					/* tp_print */
	0,					/* tp_getattr */
	0,					/* tp_setattr */
	0,					/* tp_compare */
	(reprfunc)&IdentifierWrapper::__repr__,					/* tp_repr */
	0,					/* tp_as_number */
	0,					/* tp_as_sequence */
	0,					/* tp_as_mapping */
	(hashfunc)&IdentifierWrapper::__hash__,					/* tp_hash */
	0,					/* tp_call */
	(reprfunc)&IdentifierWrapper::__str__,					/* tp_str */
	PyObject_GenericGetAttr,		/* tp_getattro */
	0,					/* tp_setattro */
	0,					/* tp_as_buffer */
	Py_TPFLAGS_HAVE_CLASS | Py_TPFLAGS_HAVE_RICHCOMPARE,/* tp_flags */
	0,					/* tp_doc */
	0,              	/* tp_traverse */
	0,					/* tp_clear */
	(richcmpfunc)&IdentifierWrapper::__richcmp__, /* tp_richcompare */
	0,					/* tp_weaklistoffset */
	0,                  /* tp_iter */
	0,                  /* tp_iternext */
	IdentifierWrapper::__methods__,		        	/* tp_methods */
	0,					/* tp_members */
    IdentifierWrapper::__getsets__, /* tp_getset */
    &PyBaseObject_Type, /* tp_base */
    0,                  /* tp_dict */
    0,                  /* tp_descr_get */
    0,                  /* tp_descr_set */
    0,                  /* tp_dictoffset S*/
    0,                  /* tp_init */
    0,                  /* tp_alloc */
    IdentifierWrapper::__new__  /*tp_new */
};

} //namespace binding

#endif /* BINDING_IDENTIFIER_HPP */
