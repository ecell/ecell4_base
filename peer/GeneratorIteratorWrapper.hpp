#ifndef PEER_GENERATOR_ITERATOR_WRAPPER_HPP
#define PEER_GENERATOR_ITERATOR_WRAPPER_HPP

#include <Python.h>
#include "generator.hpp"
#include "peer/utils.hpp"

namespace peer { namespace util {

namespace detail {

template<typename Twrapper_, typename Tgen_ = typename Twrapper_::generator_type>
struct converter_pair
{
    typedef Twrapper_ self_type;

    struct to_python_converter
    {
        static PyObject* convert(Tgen_ const& impl)
        {
            return self_type::create(impl);
        }

        static PyTypeObject* get_pytype()
        {
            return &self_type::__class__;
        }
    };

    struct to_native_converter
    {
        static void* convertible(PyObject* pyo)
        {
            if (!PyObject_TypeCheck(pyo, &self_type::__class__))
            {
                return 0;
            }
            return pyo;
        }

        static void construct(PyObject* pyo, 
                              boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            void* storage(reinterpret_cast<
                boost::python::converter::rvalue_from_python_storage<Tgen_>* >(
                    data)->storage.bytes);
            new (storage) Tgen_(reinterpret_cast<self_type*>(pyo)->impl_);
            data->convertible = storage;
        }
    };

    static void __register_converter()
    {
        util::to_native_converter<Tgen_, to_native_converter>();
        boost::python::to_python_converter<Tgen_, to_python_converter>();
    }
};

template<typename Twrapper_, typename Tgen_>
struct converter_pair<Twrapper_, ptr_generator<Tgen_> >
{
    typedef Twrapper_ self_type;

    struct to_python_converter
    {
        static PyObject* convert(Tgen_* impl)
        {
            return impl ?
                self_type::create(ptr_generator<Tgen_>(impl)):
                boost::python::incref(Py_None);
        }

        static PyTypeObject* get_pytype()
        {
            return &self_type::__class__;
        }
    };

    struct to_native_converter
    {
        static void* convertible(PyObject* pyo)
        {
            if (!PyObject_TypeCheck(pyo, &self_type::__class__))
            {
                return 0;
            }
            return pyo;
        }

        static void construct(PyObject* pyo, 
                              boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            void* storage(reinterpret_cast<
                boost::python::converter::rvalue_from_python_storage<Tgen_*>* >(
                    data)->storage.bytes);
            new (storage) Tgen_*(reinterpret_cast<self_type*>(pyo)->impl_.ptr().get());
            data->convertible = storage;
        }
    };

    static void __register_converter()
    {
        util::to_native_converter<Tgen_*, to_native_converter>();
        boost::python::to_python_converter<Tgen_*, to_python_converter>();
    }
};

} // namespace detail

template<typename Tgen_>
class GeneratorIteratorWrapper
{
    template<template<class>class TT1_, typename T1_> friend class detail::converter_pair;
public:
    typedef Tgen_ generator_type;

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

    GeneratorIteratorWrapper(Tgen_ const& impl): impl_(impl) {}

    ~GeneratorIteratorWrapper() {}

    static PyObject* create(Tgen_ const& impl)
    {
        return reinterpret_cast<PyObject*>(new GeneratorIteratorWrapper(impl));
    } 

    static void __dealloc__(GeneratorIteratorWrapper* self)
    {
        delete self;
    }

    static PyObject* next(GeneratorIteratorWrapper* self)
    {
        using namespace boost::python;

        if (!valid(self->impl_))
        {
            PyErr_SetNone(PyExc_StopIteration);
            return NULL;
        }
        return incref(object(self->impl_()).ptr());
    }


    static PyObject* __length_hint__(GeneratorIteratorWrapper* self)
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
        __name__ = static_cast<std::string>(
            extract<std::string>(object(borrowed(mod)).attr("__name__")))
            + "." + name;
        __class__.tp_name = const_cast<char*>(__name__.c_str());
        PyType_Ready(&__class__);
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

template<typename Tgen_>
inline void GeneratorIteratorWrapper<Tgen_>::__register_class(const char* name)
{
    typedef detail::converter_pair<GeneratorIteratorWrapper> base_type;
    using namespace boost::python;
    PyTypeObject* klass(GeneratorIteratorWrapper::__class_init__(name, reinterpret_cast<PyObject*>(scope().ptr())));
    Py_INCREF(klass);
    scope().attr(name) = object(borrowed(reinterpret_cast<PyObject*>(klass)));
    base_type::__register_converter(); 
}

template<typename Tgen_>
std::string GeneratorIteratorWrapper<Tgen_>::__name__;

template<typename Tgen_>
PyMethodDef GeneratorIteratorWrapper<Tgen_>::__methods__[] = {
    { "__length_hint__", (PyCFunction)GeneratorIteratorWrapper::__length_hint__, METH_NOARGS, "" },
    { NULL, NULL }
};

template<typename Tgen_>
PyTypeObject GeneratorIteratorWrapper<Tgen_>::__class__ = {
	PyObject_HEAD_INIT(&PyType_Type)
	0,					/* ob_size */
	0,                  /* tp_name */
	sizeof(GeneratorIteratorWrapper), /* tp_basicsize */
	0,					/* tp_itemsize */
	/* methods */
	(destructor)&GeneratorIteratorWrapper::__dealloc__, /* tp_dealloc */
	0,					/* tp_print */
	0,					/* tp_getattr */
	0,					/* tp_setattr */
	0,					/* tp_compare */
	0,					/* tp_repr */
	0,					/* tp_as_number */
	0,					/* tp_as_sequence */
	0,					/* tp_as_mapping */
	0,					/* tp_hash */
	0,					/* tp_call */
	0,					/* tp_str */
	PyObject_GenericGetAttr,		/* tp_getattro */
	0,					/* tp_setattro */
	0,					/* tp_as_buffer */
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER, /* tp_flags */
	0,					/* tp_doc */
	0,              	/* tp_traverse */
	0,					/* tp_clear */
	0,                  /* tp_richcompare */
	0,					/* tp_weaklistoffset */
	PyObject_SelfIter,  /* tp_iter */
	(iternextfunc)&GeneratorIteratorWrapper::next,  /* tp_iternext */
    GeneratorIteratorWrapper::__methods__   /* tp_methods */
};

} } // namespace peer::util
#endif /* PEER_GENERATOR_ITERATOR_WRAPPER_HPP */
