#ifndef PEER_EXCEPTION_HPP
#define PEER_EXCEPTION_HPP

#include <stdexcept>
#include <string>
#include <boost/python.hpp>
#include <Python.h>
#include <pyerrors.h>

namespace peer { namespace util {

template<PyObject** Vpytype_object_>
struct PyExcTraits
{
    struct unsupported_python_exception_type;
    enum { DUMMY = sizeof(unsupported_python_exception_type) };
};

#define SPECIALIZE_PYEXC_TRAITS(PTR, TYPE) \
template<> \
struct PyExcTraits<&PyExc_##PTR> \
{ \
    typedef TYPE type; \
    static PyTypeObject* pytype_object; \
}; \
\
PyTypeObject* PyExcTraits<&PyExc_##PTR>::pytype_object = reinterpret_cast<PyTypeObject*>(PyExc_##PTR);

#ifdef HAVE_PYBASEEXCEPTIONOBJECT
#define PYEXC_DICT_MEMBER_NAME dict
SPECIALIZE_PYEXC_TRAITS(BaseException, PyBaseExceptionObject)
SPECIALIZE_PYEXC_TRAITS(Exception, PyBaseExceptionObject)
SPECIALIZE_PYEXC_TRAITS(StandardError, PyBaseExceptionObject)
SPECIALIZE_PYEXC_TRAITS(LookupError, PyBaseExceptionObject)
SPECIALIZE_PYEXC_TRAITS(RuntimeError, PyBaseExceptionObject)
#else
#define PYEXC_DICT_MEMBER_NAME in_dict
SPECIALIZE_PYEXC_TRAITS(Exception, PyInstanceObject)
SPECIALIZE_PYEXC_TRAITS(StandardError, PyInstanceObject)
SPECIALIZE_PYEXC_TRAITS(LookupError, PyInstanceObject)
SPECIALIZE_PYEXC_TRAITS(RuntimeError, PyInstanceObject)
#endif

#undef SPECIALIZE_PYEXC_TRAITS

template<typename Texc_, typename TbaseTraits_ = PyExcTraits<&PyExc_StandardError> >
class ExceptionWrapper: TbaseTraits_::type
{
public:
    typedef Texc_ wrapped_type;
    typedef typename TbaseTraits_::type base_type;
    typedef TbaseTraits_ base_traits;
public:
    void* operator new(size_t)
    {
        return (*base_traits::pytype_object->tp_new)(&__class__, NULL, NULL);
    }

    void operator delete(void* ptr)
    {
        _PyObject_GC_UNTRACK(ptr);
        reinterpret_cast<PyObject*>(ptr)->ob_type->tp_free(reinterpret_cast<PyObject*>(ptr));
    }

    static PyObject* create(wrapped_type const& exc)
    {
        return reinterpret_cast<PyObject*>(new ExceptionWrapper(exc));
    }

    ExceptionWrapper(Texc_ const& impl): impl_(impl)
    {
        boost::python::tuple t(boost::python::make_tuple(impl_.what()));
        if ((*base_traits::pytype_object->tp_init)(
                reinterpret_cast<PyObject*>(this), t.ptr(), NULL) == -1)
        {
            throw std::runtime_error("Failed to initialize the base class");
        }

        this->message = boost::python::incref(
            static_cast<boost::python::object>(t[0]).ptr());
    }

    ~ExceptionWrapper()
    {
        base_traits::pytype_object->tp_clear(reinterpret_cast<PyObject*>(this));
    }

    static PyObject* __get_message__(ExceptionWrapper* self)
    {
        return self->message;
    }

    static void translate_exception(wrapped_type const& type)
    {
        PyErr_SetObject(reinterpret_cast<PyObject*>(&ExceptionWrapper::__class__), create(type));
    }

    static PyTypeObject* __class_init__(const char* name, PyObject* mod)
    {
        using namespace boost::python;
        __name__ = static_cast<std::string>(
            extract<std::string>(object(borrowed(mod)).attr("__name__")))
            + "." + name;
        __class__.tp_name = const_cast< char* >( __name__.c_str() );
        PyType_Ready(&__class__);
        return &__class__;
    }

    static void __dealloc__(ExceptionWrapper* self)
    {
        delete self;
    }

    static PyObject* __new__(PyTypeObject* type, PyObject* arg, PyObject* kwarg)
    {
        PyErr_SetString(PyExc_RuntimeError, "This class cannot be instantiated from within a Python script");
        return NULL;
    }

    static int __traverse__(ExceptionWrapper* self, visitproc visit, void *arg)
    {
        return (*base_traits::pytype_object->tp_traverse)(
            reinterpret_cast<PyObject*>(self), visit, arg);
    }

    static void __register_class(const char* name);

public:
    static PyTypeObject __class__;
    static PyGetSetDef __getsets__[];
    static std::string __name__;
    Texc_ impl_;
#ifndef HAVE_PYBASEEXCEPTIONOBJECT
    PyObject* message;
#endif 
};


template<typename Texc_, typename TbaseTraits_>
inline void ExceptionWrapper<Texc_, TbaseTraits_>::__register_class(const char* name)
{
    using namespace boost::python;
    PyTypeObject* klass(ExceptionWrapper::__class_init__(name, reinterpret_cast<PyObject*>(scope().ptr())));
    Py_INCREF(klass);
    scope().attr(name) = object(borrowed(reinterpret_cast<PyObject*>(klass)));
    register_exception_translator<Texc_>(&ExceptionWrapper::translate_exception);
}

template<typename Texc_, typename TbaseTraits_>
std::string ExceptionWrapper<Texc_, TbaseTraits_>::__name__;

template<typename Texc_, typename TbaseTraits_>
PyGetSetDef ExceptionWrapper<Texc_, TbaseTraits_>::__getsets__[] = {
    { const_cast<char*>("message"), (getter)&ExceptionWrapper::__get_message__, NULL },
    { NULL }
};

template<typename Texc_, typename TbaseTraits_>
PyTypeObject ExceptionWrapper<Texc_, TbaseTraits_>::__class__ = {
	PyObject_HEAD_INIT(NULL)
	0,					/* ob_size */
	0,                  /* tp_name */
	sizeof(ExceptionWrapper), /* tp_basicsize */
	0,					/* tp_itemsize */
	/* methods */
	(destructor)&ExceptionWrapper::__dealloc__, /* tp_dealloc */
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
	Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC | Py_TPFLAGS_BASETYPE, /* tp_flags */
	0,					/* tp_doc */
	(traverseproc)&ExceptionWrapper::__traverse__,              	/* tp_traverse */
	0,					/* tp_clear */
	0,                  /* tp_richcompare */
	0,					/* tp_weaklistoffset */
	0,                  /* tp_iter */
	0,                  /* tp_iternext */
    0,                  /* tp_methods */
    0,                  /* tp_members */
    ExceptionWrapper::__getsets__,  /* tp_getset */
    reinterpret_cast<PyTypeObject*>(base_traits::pytype_object), /* tp_base */
    0,                  /* tp_dict */
    0,                  /* tp_descr_get */
    0,                  /* tp_descr_set */
    offsetof(typename ExceptionWrapper::base_type, PYEXC_DICT_MEMBER_NAME) /* tp_dictoffset */,
    0,                  /* tp_init */
    0,                  /* tp_alloc */
    &ExceptionWrapper::__new__    /* tp_new */
};

} }

#endif /* PEER_EXCEPTION_HPP */
