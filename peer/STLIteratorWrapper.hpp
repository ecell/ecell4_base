#ifndef PEER_STLITERATORWRAPPER_HPP
#define PEER_STLITERATORWRAPPER_HPP

#include <boost/python.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>

namespace peer {

template<typename Titer_, typename Tholder_ = void*, typename Trcg_ = boost::python::return_by_value>
class STLIteratorWrapper
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
        return PyObject_New(STLIteratorWrapper, &__class__);
    }

    void operator delete(void* ptr)
    {
        reinterpret_cast<PyObject*>(ptr)->ob_type->tp_free(reinterpret_cast< PyObject*>(ptr));
    }

    template<typename Trange>
    STLIteratorWrapper(Trange const& range, Tholder_ holder = Tholder_())
        : i_(boost::begin(range)), end_(boost::end(range)), holder_(holder) {}

    template<typename Titer>
    STLIteratorWrapper(Titer const& begin, Titer const& end, Tholder_ holder = Tholder_())
        : i_(begin), end_(end), holder_(holder) {}

    ~STLIteratorWrapper()
    {
    }

public:
    static PyTypeObject* __class_init__(const char* name, PyObject* mod = 0)
    {
        if (__name__.empty())
        {
            using namespace boost::python;
            __name__ = (mod && PyModule_Check(mod) ?
                extract<std::string>(object(borrowed(mod)).attr("__name__"))()
                + ".": std::string()) + name;
            __class__.tp_name = const_cast<char*>(__name__.c_str());
            PyType_Ready(&__class__);
        }
        return &__class__;
    }

    static void __register_class(char const* name)
    {
        using namespace boost::python;
        PyTypeObject* klass(STLIteratorWrapper::__class_init__(name, reinterpret_cast<PyObject*>(scope().ptr())));
        Py_INCREF(klass);
        scope().attr(name) = object(borrowed(reinterpret_cast<PyObject*>(klass)));
    }

    template<typename Trange>
    static PyObject* create(Trange const& range, Tholder_ holder = Tholder_())
    {
        return reinterpret_cast<PyObject*>(new STLIteratorWrapper(range, holder));
    }

    static void __dealloc__(STLIteratorWrapper* self)
    {
        delete self;
    }

    static PyObject* __next__(STLIteratorWrapper* self)
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
std::string STLIteratorWrapper<Titer_, Tholder_, Trcg_>::__name__;

template<typename Titer_, typename Tholder_, typename Trcg_>
PyTypeObject STLIteratorWrapper<Titer_, Tholder_, Trcg_>::__class__ = {
    PyObject_HEAD_INIT(&PyType_Type)
    0,                    /* ob_size */
    0,                    /* tp_name */
    sizeof(STLIteratorWrapper), /* tp_basicsize */
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
    Py_TPFLAGS_HAVE_CLASS | Py_TPFLAGS_HAVE_ITER,/* tp_flags */
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

template<typename Trange_, typename Tholder_ = void*, typename Trcg_ = boost::python::return_by_value>
struct stl_iterator_range_converter
{
    typedef Trange_ native_type;

    static PyObject* convert(native_type const& v)
    {
        return reinterpret_cast<PyObject*>(STLIteratorWrapper<typename boost::range_const_iterator<native_type>::type, Tholder_, Trcg_>::create(v));
    }
};

} // namespace detail

template<typename Trange, typename Tholder, typename Trcg>
inline void register_stl_iterator_range_converter()
{
    boost::python::to_python_converter<Trange, detail::stl_iterator_range_converter<Trange, Tholder, Trcg> >();
}

template<typename Trange, typename Tholder, typename Trcg>
inline PyObject*
make_stl_iterator_wrapper(Trange const& range, Tholder holder = Tholder(), Trcg const& rcg = boost::python::return_by_value())
{
    typedef STLIteratorWrapper<typename boost::range_const_iterator<Trange>::type, Tholder, Trcg> wrapper_type;
    wrapper_type::__class_init__(typeid(wrapper_type).name(), boost::python::scope().ptr());
    return wrapper_type::create(range, holder);
}

} // namespace util

} // namespace peer

#endif /* PEER_STLITERATORWRAPPER_HPP */
