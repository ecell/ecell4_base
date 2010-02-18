#ifndef PEER_STLITERATORWRAPPER_HPP
#define PEER_STLITERATORWRAPPER_HPP

#include <boost/python.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>

namespace peer {

template< typename Titer_ >
class STLIteratorWrapper
{
protected:
    PyObject_VAR_HEAD
    Titer_ i_;
    Titer_ end_; 

public:
    static PyTypeObject __class__;

public:
    void* operator new(size_t)
    {
        return PyObject_New(STLIteratorWrapper, &__class__);
    }

    void operator delete(void* ptr)
    {
        reinterpret_cast<PyObject*>(ptr)->ob_type->tp_free(reinterpret_cast< PyObject*>(ptr));
    }

    template< typename Trange_ >
    STLIteratorWrapper(Trange_ const& range)
        : i_(boost::begin(range )), end_(boost::end(range)) {}

    STLIteratorWrapper(Titer_ const& begin, Titer_ const& end)
        : i_(begin), end_(end) {}

    ~STLIteratorWrapper()
    {
    }

public:
    static PyTypeObject* __class_init__()
    {
        PyType_Ready(&__class__);
        return &__class__;
    }

    template<typename Trange_>
    static STLIteratorWrapper* create(Trange_ const& range)
    {
        return new STLIteratorWrapper(range);
    }

    static void __dealloc__(STLIteratorWrapper* self)
    {
        delete self;
    }

    static PyObject* __next__(STLIteratorWrapper* self)
    {
        if ( self->i_ == self->end_ )
            return NULL;

        return boost::python::incref(boost::python::object(*self->i_ ++).ptr());
    }
};

template< typename Titer_ >
PyTypeObject STLIteratorWrapper< Titer_ >::__class__ = {
    PyObject_HEAD_INIT( &PyType_Type )
    0,                    /* ob_size */
    "ecell._ecs.STLIteratorWrapper", /* tp_name */
    sizeof( STLIteratorWrapper ), /* tp_basicsize */
    0,                    /* tp_itemsize */
    /* methods */
    (destructor)&STLIteratorWrapper::__dealloc__, /* tp_dealloc */
    0,                    /* tp_print */
    0,                    /* tp_getattr */
    0,                    /* tp_setattr */
    0,                    /* tp_compare */
    0,                    /* tp_repr */
    0,                    /* tp_as_number */
    0,                    /* tp_as_sequence */
    0,                    /* tp_as_mapping */
    0,                    /* tp_hash */
    0,                    /* tp_call */
    0,                    /* tp_str */
    PyObject_GenericGetAttr,        /* tp_getattro */
    0,                    /* tp_setattro */
    0,                    /* tp_as_buffer */
    Py_TPFLAGS_HAVE_CLASS | Py_TPFLAGS_HAVE_WEAKREFS | Py_TPFLAGS_HAVE_ITER,/* tp_flags */
    0,                    /* tp_doc */
    0,                    /* tp_traverse */
    0,                    /* tp_clear */
    0,                    /* tp_richcompare */
    0,                    /* tp_weaklistoffset */
    PyObject_SelfIter,  /* tp_iter */
    (iternextfunc)&STLIteratorWrapper::__next__,        /* tp_iternext */
    0,                    /* tp_methods */
    0,                    /* tp_members */
    0,                    /* tp_getset */
    0,                    /* tp_base */
    0,                    /* tp_dict */
    0,                    /* tp_descr_get */
    0,                    /* tp_descr_set */
    0,                    /* tp_dictoffset */
    0,                    /* tp_init */
    PyType_GenericAlloc,  /* tp_alloc */
    PyType_GenericNew,    /* tp_new */
    PyObject_Del,         /* tp_free */
};

namespace util {

namespace detail {

template<typename Trange_>
struct stl_iterator_range_converter
{
    typedef Trange_ native_type;

    static PyObject* convert(native_type const& v)
    {
        return reinterpret_cast<PyObject*>(STLIteratorWrapper<typename native_type::const_iterator>::create(v));
    }
};

} // namespace detail

template<typename Trange_>
inline void register_stl_iterator_range_converter()
{
    boost::python::to_python_converter<Trange_, detail::stl_iterator_range_converter<Trange_> >();
}

} // namespace util

} // namespace peer

#endif /* PEER_STLITERATORWRAPPER_HPP */
