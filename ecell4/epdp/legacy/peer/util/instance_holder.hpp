#ifndef PEER_UTIL_INSTANCE_HOLDER_HPP
#define PEER_UTIL_INSTANCE_HOLDER_HPP

#include <Python.h>
#include <typeinfo>
#include <boost/format.hpp>
#include <boost/python/handle.hpp>

namespace peer { namespace util {

template<typename Timpl_>
struct instance_holder
{
    typedef Timpl_ impl_type;

    void* operator new(size_t)
    {
        return PyObject_New(PyObject, &__class__);
    }

    void operator delete(void* ptr)
    {
        reinterpret_cast<PyObject*>(ptr)->ob_type->tp_free(reinterpret_cast<PyObject*>(ptr));
    }

    Timpl_& operator()()
    {
        return impl_;
    }

    Timpl_ const& operator()() const
    {
        return impl_;
    }

    ~instance_holder()
    {
        reinterpret_cast<impl_type*>(impl_)->~impl_type();
    }

    template<typename Tin_place_factory>
    instance_holder(Tin_place_factory const& inpf)
    {
        inpf.template apply<impl_type>(impl_);
    }

    static void __dealloc__(instance_holder* self)
    {
        delete self;
    }

    template<typename Tin_place_factory>
    static PyObject* create(Tin_place_factory const& inpf)
    {
        return reinterpret_cast<PyObject*>(new instance_holder<impl_type>(inpf));
    }

    static PyTypeObject* __class_init__(PyObject* mod = 0)
    {
        using namespace boost::python;
        if (__name__.empty())
        {
            std::string const name(
                std::string("instance_holder<")
                + typeid(impl_type).name() + ">");

            __name__ = mod && PyModule_Check(mod) ?
                extract<std::string>(object(borrowed(mod)).attr("__name__"))()
                + "." + name: name;
            __class__.tp_name = const_cast<char*>(__name__.c_str());
            PyType_Ready(&__class__);
        }
        return &__class__;
    }

protected:
    PyObject_VAR_HEAD 
    static PyTypeObject __class__;
    static std::string __name__;
    char impl_[sizeof(impl_type)];
};

template<typename Timpl_>
std::string instance_holder<Timpl_>::__name__;

template<typename Timpl_>
PyTypeObject instance_holder<Timpl_>::__class__ = {
    PyObject_HEAD_INIT(&PyType_Type)
    0,                  /* ob_size */
    0,                  /* tp_name */
    sizeof(instance_holder), /* tp_basicsize */
    0,                  /* tp_itemsize */
    /* methods */
    (destructor)&instance_holder::__dealloc__, /* tp_dealloc */
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
    Py_TPFLAGS_HAVE_CLASS | Py_TPFLAGS_BASETYPE, /* tp_flags */
    0,                  /* tp_doc */
    0,                  /* tp_traverse */
    0,                  /* tp_clear */
    0,                  /* tp_richcompare */
    0,                  /* tp_weaklistoffset */
    0,                  /* tp_iter */
    0,                  /* tp_iternext */
    0,                  /* tp_methods */
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
    PyObject_Del        /* tp_free */
};

template<typename Timpl, typename Tin_place_factory>
inline void install_instance_holder(PyObject* obj, Tin_place_factory const& inpf)
{
    instance_holder<Timpl>::__class_init__();

    boost::python::handle<> holder(instance_holder<Timpl>::create(inpf));
    if (PyObject_SetAttrString(obj, "__instance_holder__", holder.get()))
    {
        PyErr_Clear();
        boost::python::handle<> repr(PyObject_Repr(obj));
        throw std::invalid_argument(
            (boost::format("object %s is not assignable") %
                PyString_AS_STRING(repr.get())).str());
    }
}

} } // namespace peer::util
#endif /* PEER_UTIL_INSTANCE_HOLDER_HPP */
