#ifndef PEER_WRAPPERS_ITERATOR_STL_ITERATOR_WRAPPER_HPP
#define PEER_WRAPPERS_ITERATOR_STL_ITERATOR_WRAPPER_HPP

#include <boost/python.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>

namespace peer { namespace wrappers {

template<typename Titer_, typename Tholder_ = void*, typename Trcg_ = boost::python::return_by_value>
class stl_iterator_wrapper
{
    typedef typename boost::iterator_reference<Titer_>::type iterator_reference;
    typedef typename Trcg_::template apply<iterator_reference>::type result_converter_type;

protected:
    PyObject_VAR_HEAD
    Titer_ i_;
    Titer_ end_;
    Tholder_ holder_;

public:
    static PyTypeObject __class__;
    static std::string __name__;

public:
    void* operator new(size_t)
    {
        return PyObject_New(stl_iterator_wrapper, &__class__);
    }

    void operator delete(void* ptr)
    {
        reinterpret_cast<PyObject*>(ptr)->ob_type->tp_free(reinterpret_cast< PyObject*>(ptr));
    }

    template<typename Trange>
    stl_iterator_wrapper(Trange const& range, Tholder_ holder = Tholder_())
        : i_(boost::begin(range)), end_(boost::end(range)), holder_(holder) {}

    template<typename Titer>
    stl_iterator_wrapper(Titer const& begin, Titer const& end, Tholder_ holder = Tholder_())
        : i_(begin), end_(end), holder_(holder) {}

    ~stl_iterator_wrapper()
    {
    }

public:
    static PyTypeObject* __class_init__(const char* name, PyObject* mod = 0)
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

    static void __register_class(char const* name)
    {
        using namespace boost::python;
        PyTypeObject* klass(stl_iterator_wrapper::__class_init__(name, reinterpret_cast<PyObject*>(scope().ptr())));
        Py_INCREF(klass);
        scope().attr(name) = object(borrowed(reinterpret_cast<PyObject*>(klass)));
    }

    template<typename Trange>
    static PyObject* create(Trange const& range, Tholder_ holder = Tholder_())
    {
        return reinterpret_cast<PyObject*>(new stl_iterator_wrapper(range, holder));
    }

    static void __dealloc__(stl_iterator_wrapper* self)
    {
        delete self;
    }

    static PyObject* __next__(stl_iterator_wrapper* self)
    {
        if (self->i_ == self->end_)
            return NULL;

        try
        {
            return result_converter_type()(*self->i_ ++);
        }
        catch (boost::python::error_already_set const&)
        {
            return NULL;
        }
    }
};

template<typename Titer_, typename Tholder_, typename Trcg_>
std::string stl_iterator_wrapper<Titer_, Tholder_, Trcg_>::__name__;

template<typename Titer_, typename Tholder_, typename Trcg_>
PyTypeObject stl_iterator_wrapper<Titer_, Tholder_, Trcg_>::__class__ = {
    PyObject_HEAD_INIT(&PyType_Type)
    0,                    /* ob_size */
    0,                    /* tp_name */
    sizeof(stl_iterator_wrapper), /* tp_basicsize */
    0,                    /* tp_itemsize */
    /* methods */
    (destructor)&stl_iterator_wrapper::__dealloc__, /* tp_dealloc */
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
    Py_TPFLAGS_HAVE_CLASS | Py_TPFLAGS_HAVE_ITER,/* tp_flags */
    0,                    /* tp_doc */
    0,                    /* tp_traverse */
    0,                    /* tp_clear */
    0,                    /* tp_richcompare */
    0,                    /* tp_weaklistoffset */
    PyObject_SelfIter,  /* tp_iter */
    (iternextfunc)&stl_iterator_wrapper::__next__,        /* tp_iternext */
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

template<typename Trange, typename Tholder, typename Trcg>
inline PyObject*
make_stl_iterator_wrapper(Trange const& range, Tholder holder = Tholder(), Trcg const& rcg = boost::python::return_by_value())
{
    typedef stl_iterator_wrapper<typename boost::range_const_iterator<Trange>::type, Tholder, Trcg> wrapper_type;
    wrapper_type::__class_init__(typeid(wrapper_type).name(), boost::python::scope().ptr());
    return wrapper_type::create(range, holder);
}

} } // namespace peer::wrappers

#endif /* PEER_WRAPPERS_ITERATOR_STL_ITERATOR_WRAPPER_HPP */
