#ifndef PEER_STLCONTAINERWRAPPER_HPP
#define PEER_STLCONTAINERWRAPPER_HPP

#include <boost/python.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size_type.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/reference.hpp>
#include <boost/range/const_iterator.hpp>

#include "peer/STLIteratorWrapper.hpp"

namespace peer {

namespace detail {

template<typename T_>
struct default_container_wrapper_policy
{
    typedef typename boost::range_size<T_>::type size_type;
    typedef typename boost::range_value<T_>::type value_type;
    typedef value_type& reference;
    typedef value_type const& const_reference;
    typedef typename boost::range_iterator<T_>::type iterator;
    typedef typename boost::range_const_iterator<T_>::type const_iterator;

    static size_type size(T_ const& c)
    {
        return boost::size(c);
    }

    static void set(T_& c, size_type i, const_reference v)
    {
        c[i] = v;
    }

    static const_reference get(T_ const& c, size_type i)
    {
        return c[i];
    }

    static iterator begin(T_& c)
    {
        return boost::begin(c);
    }

    static const_iterator begin(T_ const& c)
    {
        return boost::begin(c);
    }

    static iterator end(T_& c)
    {
        return boost::end(c);
    }

    static const_iterator end(T_ const& c)
    {
        return boost::end(c);
    }
};

} // namespace detail

template<typename Timpl_, typename Tholder_, template<typename> class TTpolicy_ = detail::default_container_wrapper_policy>
class STLContainerWrapper
{
private:
    typedef Timpl_ impl_type;
    typedef TTpolicy_<impl_type> policy_type;
    typedef STLIteratorWrapper<
            typename policy_type::const_iterator,
            boost::python::handle<> > iterator_wrapper_type;
protected:
    PyObject_VAR_HEAD
    Tholder_ impl_;

public:
    static PyTypeObject __class__;
    static PySequenceMethods __sequence_methods__;
    static std::string __name__;

public:
    void* operator new(size_t)
    {
        return PyObject_New(STLContainerWrapper, &__class__);
    }

    void operator delete(void* ptr)
    {
        reinterpret_cast<PyObject*>(ptr)->ob_type->tp_free(reinterpret_cast< PyObject*>(ptr));
    }

    // STLContainerWrapper(Tholder_ const& impl): impl_(impl) {}

    STLContainerWrapper(Tholder_ impl): impl_(impl) {}

    ~STLContainerWrapper()
    {
    }

public:
    static Py_ssize_t __sq_len__(STLContainerWrapper const* self)
    {
        return policy_type::size(*self->impl_);
    }

    static PyObject* __sq_item__(STLContainerWrapper const* self, Py_ssize_t idx)
    {
        if (idx < 0 || idx >= static_cast<Py_ssize_t>(policy_type::size(*self->impl_)))
        {
            PyErr_Format(PyExc_IndexError, "index out of range: %zd", idx);
            return NULL;
        }
        try
        {
            return boost::python::incref(boost::python::object(policy_type::get(*self->impl_, idx)).ptr());
        }
        catch (boost::python::error_already_set const&)
        {
        }
        return NULL;
    }

    static int __sq_ass_item__(STLContainerWrapper* self, Py_ssize_t idx, PyObject *val)
    {
        if (idx < 0 || idx >= static_cast<Py_ssize_t>(policy_type::size(*self->impl_)))
        {
            PyErr_Format(PyExc_IndexError, "index out of range: %zd", idx);
            return -1;
        }

        try
        {
            policy_type::set(*self->impl_, idx, boost::python::extract<typename policy_type::value_type>(val)());
        }
        catch (boost::python::error_already_set const&)
        {
            return NULL;
        }
        return 0;
    }

    static int __sq_contains__(STLContainerWrapper const* self, PyObject *val)
    {
        boost::python::extract<typename policy_type::value_type> _val(val);
        if (!_val.check())
        {
            return 0;
        }

        typename policy_type::const_iterator e(policy_type::end(static_cast<impl_type const&>(*self->impl_)));
        return e != std::find(policy_type::begin(static_cast<impl_type const&>(*self->impl_)), e, _val());
    }

    static PyObject* __iter__(STLContainerWrapper const* self)
    {
        return iterator_wrapper_type::create(*self->impl_,
            boost::python::handle<>(boost::python::borrowed(
                const_cast<PyObject*>(
                    reinterpret_cast<PyObject const*>(self)))));
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

    static void __register_class(char const *name)
    {
        using namespace boost::python;
        PyObject* mod(scope().ptr());
        iterator_wrapper_type::__class_init__(
                (std::string(name) + ".Iterator").c_str(), mod);
        PyTypeObject* klass(STLContainerWrapper::__class_init__(name, mod));
        Py_INCREF(klass);
        scope().attr(name) = object(borrowed(reinterpret_cast<PyObject*>(klass)));
    }

    static STLContainerWrapper* create(Tholder_ impl)
    {
        return new STLContainerWrapper(impl);
    }

    static void __dealloc__(STLContainerWrapper* self)
    {
        delete self;
    }
};

template<typename Timpl_, typename Tholder_, template<typename>class TTpolicy_>
PySequenceMethods STLContainerWrapper<Timpl_, Tholder_, TTpolicy_>::__sequence_methods__ = {
    (lenfunc) &STLContainerWrapper::__sq_len__,         /* sq_length */
    (binaryfunc) 0,                                     /* sq_concat */
    (ssizeargfunc) 0,                                   /* sq_repeat */
    (ssizeargfunc) &STLContainerWrapper::__sq_item__,   /* sq_item */
    (ssizessizeargfunc) 0,                              /* sq_slice */
    (ssizeobjargproc) &STLContainerWrapper::__sq_ass_item__,    /* sq_ass_item */
    (ssizessizeobjargproc) 0,                           /* sq_ass_slice */
    (objobjproc) &STLContainerWrapper::__sq_contains__, /* sq_contains */
    (binaryfunc) 0,                                     /* sq_inplace_concat */
    (ssizeargfunc) 0,                                   /* sq_inplace_repeat */
};

template<typename Titer_, typename Tholder_, template<typename>class TTpolicy_>
std::string STLContainerWrapper<Titer_, Tholder_, TTpolicy_>::__name__;

template<typename Titer_, typename Tholder_, template<typename>class TTpolicy_>
PyTypeObject STLContainerWrapper<Titer_, Tholder_, TTpolicy_>::__class__ = {
    PyObject_HEAD_INIT(&PyType_Type)
    0,                    /* ob_size */
    0,                    /* tp_name */
    sizeof(STLContainerWrapper), /* tp_basicsize */
    0,                    /* tp_itemsize */
    /* methods */
    (destructor)&STLContainerWrapper::__dealloc__, /* tp_dealloc */
    0,                    /* tp_print */
    0,                    /* tp_getattr */
    0,                    /* tp_setattr */
    0,                    /* tp_compare */
    0,                    /* tp_repr */
    0,                    /* tp_as_number */
    &STLContainerWrapper::__sequence_methods__,  /* tp_as_sequence */
    0,                    /* tp_as_mapping */
    0,                    /* tp_hash */
    0,                    /* tp_call */
    0,                    /* tp_str */
    PyObject_GenericGetAttr,        /* tp_getattro */
    0,                    /* tp_setattro */
    0,                    /* tp_as_buffer */
    Py_TPFLAGS_HAVE_CLASS | Py_TPFLAGS_HAVE_ITER, /* tp_flags */
    0,                    /* tp_doc */
    0,                    /* tp_traverse */
    0,                    /* tp_clear */
    0,                    /* tp_richcompare */
    0,                    /* tp_weaklistoffset */
    (getiterfunc)&STLContainerWrapper::__iter__, /* tp_iter */
    0,                    /* tp_iternext */
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

} // namespace peer

#endif /* PEER_STLCONTAINERWRAPPER_HPP */
