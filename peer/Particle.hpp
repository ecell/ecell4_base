#ifndef PEER_PARTICLE_HPP
#define PEER_PARTICLE_HPP

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
#include <iostream>

#include "numpy/type_mappings.hpp"
#include "pickle_support.hpp"

namespace peer {

template<typename Timpl_>
class ParticleWrapper
{
public:
    struct to_python_converter
    {
        static PyObject* convert(Timpl_ const& impl)
        {
            return ParticleWrapper::create(impl);
        }

        static PyTypeObject* get_pytype()
        {
            return &ParticleWrapper::__class__;
        }
    };

    struct to_native_converter
    {
        static void* convertible(PyObject* pyo)
        {
            if (!PyObject_TypeCheck(pyo, &ParticleWrapper::__class__))
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
            new (storage) Timpl_(reinterpret_cast<ParticleWrapper*>(pyo)->impl_);
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
        __class__.tp_name = __name__.c_str();
        PyType_Ready(&__class__);
        dict(borrowed(__class__.tp_dict))["__safe_for_unpickling__"] = object(true);
        return &__class__;
    }

    static PyObject* create()
    {
        return reinterpret_cast<PyObject*>(new ParticleWrapper());
    }

    static PyObject* create(Timpl_ const& impl)
    {
        return reinterpret_cast<PyObject*>(new ParticleWrapper(impl));
    }

    static PyObject* __new__(PyTypeObject* klass, PyObject* arg, PyObject* kwarg)
    {
        PyObject* retval = NULL;
        switch (PyTuple_Size(arg))
        {
        default:
            PyErr_SetString(PyExc_TypeError, "the number of arguments must be either 0 or 3");
            return NULL;

        case 4:
            retval = create();
            if (set_position(reinterpret_cast<ParticleWrapper*>(retval), 
                        PyTuple_GetItem(arg, 0), 0)
                || set_radius(reinterpret_cast<ParticleWrapper*>(retval), 
                        PyTuple_GetItem(arg, 1), 0)
                || set_D(reinterpret_cast<ParticleWrapper*>(retval), 
                        PyTuple_GetItem(arg, 2), 0)
                || set_sid(reinterpret_cast<ParticleWrapper*>(retval), 
                        PyTuple_GetItem(arg, 3), 0))
            {
                ParticleWrapper::operator delete(retval);
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

    static PyObject* __str__(ParticleWrapper* self)
    {
        std::string retval(boost::lexical_cast<std::string>(self->impl_));
        return PyString_FromStringAndSize(retval.data(), retval.size());
    }

    static PyObject* __repr__(ParticleWrapper* self)
    {
        return __str__(self);
    }

    static void __dealloc__(ParticleWrapper* self)
    {
        delete self;
    }

    static PyObject* get_position(ParticleWrapper* self)
    {
        typename Timpl_::position_type const& pos(self->impl_.position());
        const npy_intp dims[1] = { boost::size(pos) };
        PyObject* retval = PyArray_New(&PyArray_Type, 1,
                const_cast<npy_intp*>(dims),
                util::get_numpy_typecode<typename Timpl_::length_type>::value,
                NULL,
                NULL,
                0, NPY_CARRAY, NULL);
        if (retval == NULL)
            return NULL;
        std::memmove(PyArray_DATA(retval), &pos[0],
                sizeof(typename Timpl_::length_type) * boost::size(pos));
        return retval;
    }

    static int set_position(ParticleWrapper* self, PyObject* val, void *)
    {
        if (!PySequence_Check(val))
        {
            PyErr_SetString(PyExc_TypeError, "argument must be a sequence");
            return -1;
        }
        
        if (PySequence_Size(val) != 3)
        {
            PyErr_SetString(PyExc_ValueError, "argument must be a sequence of 3 elements");
            return -1;
        }

        PyObject* items[3] = {
            PySequence_GetItem(val, 0),
            PySequence_GetItem(val, 1),
            PySequence_GetItem(val, 2)
        };
        const typename Timpl_::position_type tmp(
            PyFloat_AsDouble(items[0]),
            PyFloat_AsDouble(items[1]),
            PyFloat_AsDouble(items[2]));
        Py_XDECREF(items[0]);
        Py_XDECREF(items[1]);
        Py_XDECREF(items[2]);
        if (PyErr_Occurred())
            return -1;
        self->impl_.position() = tmp;
        return 0;
    }

    static PyObject* get_radius(ParticleWrapper* self)
    {
        return PyFloat_FromDouble(self->impl_.radius());
    }

    static int set_radius(ParticleWrapper* self, PyObject* val, void *)
    {
        const double tmp(PyFloat_AsDouble(val));
        if (PyErr_Occurred())
            return -1;
        self->impl_.radius() = tmp;
        return 0;
    }

    static PyObject* get_D(ParticleWrapper* self)
    {
        return PyFloat_FromDouble(self->impl_.D());
    }

    static int set_D(ParticleWrapper* self, PyObject* val, void *)
    {
        const double tmp(PyFloat_AsDouble(val));
        if (PyErr_Occurred())
            return -1;
        self->impl_.D() = tmp;
        return 0;
    }

    static PyObject* get_sid(ParticleWrapper* self)
    {
        return boost::python::incref(
            boost::python::object(self->impl_.sid()).ptr());
    }

    static int set_sid(ParticleWrapper* self, PyObject* val, void *)
    try
    {
        self->impl_.sid() = boost::python::extract<typename Timpl_::species_id_type>(val);
        return 0;
    }
    catch (boost::python::error_already_set const&)
    {
        return -1;
    }

    static PyObject* __getstate__(ParticleWrapper* self)
    try
    {
        return boost::python::incref(
            boost::python::make_tuple(
                boost::python::borrowed(get_position(self)),
                boost::python::borrowed(get_radius(self)),
                boost::python::borrowed(get_sid(self))).ptr());
    }
    catch (boost::python::error_already_set const&)
    {
        return NULL;
    }

    static PyObject* __reduce__(ParticleWrapper* self)
    {
        return pickle::reduce(reinterpret_cast<PyObject*>(self));
    }

    static PyObject* __reduce_ex__(ParticleWrapper* self, PyObject* arg)
    {
        return pickle::reduce(reinterpret_cast<PyObject*>(self));
    }

    static long __hash__(ParticleWrapper* self)
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

    static Py_ssize_t __sq_len__(PyObject *self)
    {
        return 3;
    }

    static PyObject* __sq_item__(PyObject *self, Py_ssize_t idx)
    {
        if (idx < 0 || idx >= 3)
        {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            return NULL;
        }

        switch (idx)
        {
        case 0:
            return get_position(reinterpret_cast<ParticleWrapper*>(self));
        case 1:
            return get_radius(reinterpret_cast<ParticleWrapper*>(self));
        case 2:
            return get_D(reinterpret_cast<ParticleWrapper*>(self));
        case 3:
            return get_sid(reinterpret_cast<ParticleWrapper*>(self));
        }
        return NULL; // never get here
    }

    static int __sq_ass_item__(PyObject *self, Py_ssize_t idx, PyObject *val)
    {
        if (idx < 0 || idx >= 4)
        {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            return -1;
        }

        switch (idx)
        {
        case 0:
            return set_position(reinterpret_cast<ParticleWrapper*>(self), val, 0);
        case 1:
            return set_radius(reinterpret_cast<ParticleWrapper*>(self), val, 0);
        case 2:
            return set_D(reinterpret_cast<ParticleWrapper*>(self), val, 0);
        case 3:
            return set_sid(reinterpret_cast<ParticleWrapper*>(self), val, 0);
        }

        return -1; // never get here
    }

    // this method borrows a reference to "name" from the caller.
    // i.e. name should be statically allocated
    static void __register_class(char const* name)
    {
        using namespace boost::python;
        pickle::register_reconstructor();
        PyTypeObject* klass(ParticleWrapper::__class_init__(name, reinterpret_cast<PyObject*>(scope().ptr())));
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
    ParticleWrapper(): impl_() {}

    ParticleWrapper(Timpl_ const& impl): impl_(impl) {}

    void* operator new(size_t)
    {
        PyObject* retval = PyObject_New(PyObject, &__class__);
        return retval;
    }

    void operator delete(void* ptr)
    {
        reinterpret_cast<PyObject*>(ptr)->ob_type->tp_free(reinterpret_cast<PyObject*>(ptr));
    }

    ~ParticleWrapper() {}

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
std::string ParticleWrapper<Timpl_>::__name__;

template<typename Timpl_>
PyMethodDef ParticleWrapper<Timpl_>::__methods__[] = {
    { "__getstate__", (PyCFunction)ParticleWrapper::__getstate__, METH_NOARGS, "" },
    { "__reduce__", (PyCFunction)ParticleWrapper::__reduce__, METH_NOARGS, "" },
    { "__reduce_ex__", (PyCFunction)ParticleWrapper::__reduce_ex__, METH_O, "" },
    { NULL, NULL }
};

template<typename Timpl_>
PyGetSetDef ParticleWrapper<Timpl_>::__getsets__[] = {
    {
        const_cast<char*>("position"),
        (getter)ParticleWrapper::get_position,
        (setter)ParticleWrapper::set_position,
        const_cast<char*>("")
    },
    {
        const_cast<char*>("radius"),
        (getter)ParticleWrapper::get_radius,
        (setter)ParticleWrapper::set_radius,
        const_cast<char*>("")
    },
    {
        const_cast<char*>("D"),
        (getter)ParticleWrapper::get_D,
        (setter)ParticleWrapper::set_D,
        const_cast<char*>("")
    },
    {
        const_cast<char*>("sid"),
        (getter)ParticleWrapper::get_sid,
        (setter)ParticleWrapper::set_sid,
        const_cast<char*>("")
    },
    { NULL }
};

template<typename Timpl_>
PySequenceMethods ParticleWrapper<Timpl_>::__sequence_methods__ = {
    (lenfunc) &ParticleWrapper::__sq_len__,             /* sq_length */
    (binaryfunc) 0,                                     /* sq_concat */
    (ssizeargfunc) 0,                                   /* sq_repeat */
    (ssizeargfunc) &ParticleWrapper::__sq_item__,       /* sq_item */
    (ssizessizeargfunc) 0,                              /* sq_slice */
    (ssizeobjargproc) &ParticleWrapper::__sq_ass_item__,    /* sq_ass_item */
    (ssizessizeobjargproc) 0,                           /* sq_ass_slice */
    (objobjproc) 0,                                     /* sq_contains */
    (binaryfunc) 0,                                     /* sq_inplace_concat */
    (ssizeargfunc) 0,                                   /* sq_inplace_repeat */
};

template<typename Timpl_>
PyTypeObject ParticleWrapper<Timpl_>::__class__ = {
	PyObject_HEAD_INIT(&PyType_Type)
	0,					/* ob_size */
	0,                  /* tp_name */
	sizeof(ParticleWrapper), /* tp_basicsize */
	0,					/* tp_itemsize */
	/* methods */
	(destructor)&ParticleWrapper::__dealloc__, /* tp_dealloc */
	0,					/* tp_print */
	0,					/* tp_getattr */
	0,					/* tp_setattr */
	0,					/* tp_compare */
	(reprfunc)&ParticleWrapper::__repr__,					/* tp_repr */
	0,					/* tp_as_number */
	&ParticleWrapper::__sequence_methods__,	/* tp_as_sequence */
	0,					/* tp_as_mapping */
	(hashfunc)&ParticleWrapper::__hash__,					/* tp_hash */
	0,					/* tp_call */
	(reprfunc)&ParticleWrapper::__str__,					/* tp_str */
	PyObject_GenericGetAttr,		/* tp_getattro */
	0,					/* tp_setattro */
	0,					/* tp_as_buffer */
	Py_TPFLAGS_HAVE_CLASS | Py_TPFLAGS_HAVE_RICHCOMPARE,/* tp_flags */
	0,					/* tp_doc */
	0,              	/* tp_traverse */
	0,					/* tp_clear */
	0,                  /* tp_richcompare */
	0,					/* tp_weaklistoffset */
	0,                  /* tp_iter */
	0,                  /* tp_iternext */
	ParticleWrapper::__methods__,		        	/* tp_methods */
	0,					/* tp_members */
    ParticleWrapper::__getsets__, /* tp_getset */
    &PyBaseObject_Type, /* tp_base */
    0,                  /* tp_dict */
    0,                  /* tp_descr_get */
    0,                  /* tp_descr_set */
    0,                  /* tp_dictoffset */
    0,                  /* tp_init */
    0,                  /* tp_alloc */
    ParticleWrapper::__new__,  /*tp_new */
    0                   /* tp_free */
};

} //namespace peer

#endif /* PEER_PARTICLE_HPP */
