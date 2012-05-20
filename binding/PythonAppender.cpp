#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <map>
#include <cstdio>
#include <boost/foreach.hpp>
#include "PythonAppender.hpp"
#include "binding_common.hpp"

namespace binding {

typedef std::map<enum Logger::level, boost::python::object> loglevel_map_type;

static boost::python::object logging_module;
static loglevel_map_type loglevel_map;

static void import_logging_module()
{
    if (logging_module.ptr() != Py_None)
        return;

    logging_module = boost::python::object(boost::python::borrowed(PyImport_Import(boost::python::handle<>(PyString_FromString("logging")).get())));
    if (PyErr_Occurred())
    {
        boost::python::throw_error_already_set();
    }

    loglevel_map[Logger::L_OFF]     = getattr(logging_module, "NOTSET");
    loglevel_map[Logger::L_DEBUG]   = getattr(logging_module, "DEBUG");
    loglevel_map[Logger::L_INFO]    = getattr(logging_module, "INFO");
    loglevel_map[Logger::L_WARNING] = getattr(logging_module, "WARNING");
    loglevel_map[Logger::L_ERROR]   = getattr(logging_module, "ERROR");
    loglevel_map[Logger::L_FATAL]   = getattr(logging_module, "CRITICAL");
}

class PythonAppender: public LogAppender
{
public:
    virtual ~PythonAppender() {}

    virtual void operator()(enum Logger::level lv,
                            char const* name, char const** chunks)
    {
        std::string msg;
        for (char const** p = chunks; *p; ++p)
            msg.append(*p);
        handle_(makeRecord_(name, lv, "", 0, msg.c_str(), NULL, NULL, NULL));
    }

    virtual void flush()
    {
        flush_();
    }

    PythonAppender(boost::python::object handler)
        : handler_(handler), flush_(handler_.attr("flush")),
          handle_(handler_.attr("handle")),
          makeRecord_(handler_.attr("makeRecord")) {}  

private:
    boost::python::object handler_;
    boost::python::object flush_;
    boost::python::object handle_;
    boost::python::object makeRecord_;
};

class CppLoggerHandler
{
public:
// XXX: logging.Handler is a old-style class... *sigh*
#ifdef HAVE_OLD_FASHIONED_LOGGER_CLASSES
    static PyObject* __class_init__(const char* name, PyObject* mod)
    {
        using namespace boost::python;

        import_logging_module();
        object mod_name(object(borrowed(mod)).attr("__name__"));

        boost::python::handle<> dict(PyDict_New());
        if (PyDict_SetItemString(dict.get(), "__module__", mod_name.ptr()))
            return NULL;

        boost::python::object bases(
            boost::python::make_tuple(
                getattr(logging_module, "Handler")));

        boost::python::object _name(name);

        boost::python::handle<> klass(
            PyClass_New(bases.ptr(), dict.get(), _name.ptr()));
        if (PyErr_Occurred())
            return NULL;

        for (PyMethodDef* meth = __methods__; meth->ml_name; ++meth)
        {
            boost::python::handle<> _meth(
                boost::python::allow_null(
                    meth->ml_flags & METH_STATIC ?
                        PyStaticMethod_New(PyCFunction_New(meth, NULL)):
                        PyDescr_NewMethod(&PyInstance_Type, meth)));
            if (!_meth)
                return NULL;
            if (PyDict_SetItemString(dict.get(), meth->ml_name, _meth.get()))
                return NULL;
        }

        for (PyGetSetDef* gsp = __getsets__; gsp->name; ++gsp)
        {
            boost::python::handle<> descr(
                boost::python::allow_null(
                    PyDescr_NewGetSet(&PyInstance_Type, gsp)));
            if (!descr)
                return NULL;
            if (PyDict_SetItemString(dict.get(), gsp->name, descr.get()))
                return NULL;
        }

        __class__ = klass;

        return klass.get();
    }

    static CppLoggerHandler* get_self(PyObject* _self)
    {
        boost::python::handle<> ptr(
            boost::python::allow_null(
                PyObject_GetAttrString(_self, "__ptr__")));
        if (!ptr)
            return 0;
        if (!PyCObject_Check(ptr.get()))
        {
            PyErr_Format(PyExc_TypeError,
                    "self.__ptr__ must be a PyCObject instance (got %s)",
                    ptr.get()->ob_type->tp_name);
            return 0;
        }
        return reinterpret_cast<CppLoggerHandler*>(PyCObject_AsVoidPtr(ptr.get()));
    }

    static PyObject* __init__(PyObject* _self, PyObject* logger)
    try
    {
        BOOST_ASSERT(PyInstance_Check(_self));

        CppLoggerHandler* self(new CppLoggerHandler(
            *boost::python::extract<Logger*>(logger)()));

        boost::python::handle<> ptr(
            boost::python::allow_null(
                PyCObject_FromVoidPtr(self, &__dealloc__)));
        if (!ptr)
        {
            delete self;
            return NULL;
        }

        if (PyObject_SetAttrString(_self, "__ptr__", ptr.get()))
            return NULL;

        boost::python::handle<> bases(
            boost::python::borrowed(
                reinterpret_cast<PyInstanceObject*>(_self)->in_class->cl_bases));
   
        if (!PyTuple_Check(bases.get()))
            return NULL;

        boost::python::handle<> super_arg(PyTuple_Pack(1, _self));

        for (std::size_t i(0), l(PyTuple_GET_SIZE(bases.get())); i < l; ++i)
        {
            boost::python::handle<> base(
                boost::python::borrowed(PyTuple_GET_ITEM(bases.get(), i)));
            boost::python::handle<> super_init(
                boost::python::borrowed(
                    PyObject_GetAttrString(base.get(), "__init__")));
            if (PyErr_Occurred())
            {
                PyErr_Clear();
                continue;
            }

            if (!boost::python::handle<>(
                boost::python::allow_null(
                    PyObject_Call(super_init.get(), super_arg.get(), NULL))))

            {
                return NULL; 
            }
        }

        return boost::python::incref(Py_None);
    }
    catch (boost::python::error_already_set const&)
    {
        return NULL;
    }

    static PyObject* __init__(PyObject* _self, PyObject* logger)
    try
    {
        BOOST_ASSERT(PyInstance_Check(_self));

        CppLoggerHandler* self(new CppLoggerHandler(
            *boost::python::extract<Logger*>(logger)()));

        boost::python::handle<> ptr(
            boost::python::allow_null(
                PyCObject_FromVoidPtr(self, &__dealloc__)));
        if (!ptr)
        {
            delete self;
            return NULL;
        }

        if (PyObject_SetAttrString(_self, "__ptr__", ptr.get()))
            return NULL;

        boost::python::handle<> bases(
            boost::python::borrowed(
                reinterpret_cast<PyInstanceObject*>(_self)->in_class->cl_bases));
   
        if (!PyTuple_Check(bases.get()))
            return NULL;

        boost::python::handle<> super_arg(PyTuple_Pack(1, _self));

        for (std::size_t i(0), l(PyTuple_GET_SIZE(bases.get())); i < l; ++i)
        {
            boost::python::handle<> base(
                boost::python::borrowed(PyTuple_GET_ITEM(bases.get(), i)));
            boost::python::handle<> super_init(
                boost::python::borrowed(
                    PyObject_GetAttrString(base.get(), "__init__")));
            if (PyErr_Occurred())
            {
                PyErr_Clear();
                continue;
            }

            if (!boost::python::handle<>(
                boost::python::allow_null(
                    PyObject_Call(super_init.get(), super_arg.get(), NULL))))

            {
                return NULL; 
            }
        }

        return boost::python::incref(Py_None);
    }
    catch (boost::python::error_already_set const&)
    {
        return NULL;
    }
#else
    void* operator new(size_t)
    {
        PyObject* retval = PyObject_New(PyObject, &__class__);
        return retval;
    }

    void operator delete(void* ptr)
    {
        reinterpret_cast<PyObject*>(ptr)->ob_type->tp_free(reinterpret_cast<PyObject*>(ptr));
    }

    static PyObject* __class_init__(const char* name, PyObject* mod)
    {
        using namespace boost::python;

        import_logging_module();
        __name__ = static_cast<std::string>(
            extract<std::string>(object(borrowed(mod)).attr("__name__")))
            + "." + name;
        __class__.tp_name = const_cast<char*>(__name__.c_str());
        __class__.tp_base = reinterpret_cast<PyTypeObject*>(getattr(logging_module, "Handler").ptr());
        if (PyType_Ready(&__class__) < 0)
        {
            return 0;
        }
        return reinterpret_cast<PyObject*>(&__class__);
    }

    static CppLoggerHandler* get_self(PyObject* _self)
    {
        return reinterpret_cast<CppLoggerHandler*>(_self);
    }

    static PyObject* __new__(PyTypeObject* klass, PyObject* arg, PyObject* kwarg)
    try
    {
        if (PyTuple_Size(arg) != 1)
        {
            PyErr_SetString(PyExc_TypeError, "the number of arguments must be 1");
            return NULL;
        }

        CppLoggerHandler* retval(
            new CppLoggerHandler(
                *boost::python::extract<Logger*>(
                    PyTuple_GetItem(arg, 0))()));

        if (PyErr_Occurred())
        {
            boost::python::decref(reinterpret_cast<PyObject*>(retval));
            return NULL;
        }
        return reinterpret_cast<PyObject*>(retval);
    }
    catch (boost::python::error_already_set const&)
    {
        return NULL;
    }

    static int __init__(PyObject* self, PyObject* arg, PyObject* kwarg)
    {
        using namespace boost::python;
        return handle<>(allow_null(
            PyObject_Call(
                PyObject_GetAttrString(
                    reinterpret_cast<PyObject*>(Py_TYPE(self)->tp_base),
                    const_cast<char*>("__init__")),
                boost::python::handle<>(PyTuple_Pack(1, self)).get(),
                NULL))) ? 0: 1;
    }

    static int __traverse__(PyObject *self, visitproc visit, void *arg)
    {
        Py_VISIT(reinterpret_cast<CppLoggerHandler*>(self)->__dict__.get());
        return 0;
    }

#endif
    ~CppLoggerHandler()
    {
        if (__weakreflist__)
        {
            PyObject_ClearWeakRefs(reinterpret_cast<PyObject*>(this));
        }
    }

    static void __dealloc__(void* ptr)
    {
        delete reinterpret_cast<CppLoggerHandler*>(ptr);
    }

    static PyObject* createLock(PyObject* _self)
    {
        CppLoggerHandler* const self(get_self(_self));
        if (!self)
            return NULL;
        return boost::python::incref(Py_None);
    }

    static PyObject* acquire(PyObject* _self)
    {
        CppLoggerHandler* const self(get_self(_self));
        if (!self)
            return NULL;
        return boost::python::incref(Py_None);
    }

    static PyObject* release(PyObject* _self)
    {
        CppLoggerHandler* const self(get_self(_self));
        if (!self)
            return NULL;
        return boost::python::incref(Py_None);
    }

    static PyObject* setLevel(PyObject* _self, PyObject* _level)
    {
        CppLoggerHandler* const self(get_self(_self));
        if (!self)
            return NULL;
        enum Logger::level const level(translate_level(_level));
        if (level == Logger::L_OFF)
        {
            PyErr_SetString(PyExc_ValueError, "invalid loglevel");
            return NULL;
        }
        self->impl_.level(level);
        return boost::python::incref(Py_None);
    }

    static PyObject* getLevel(PyObject* _self)
    {
        CppLoggerHandler* const self(get_self(_self));
        if (!self)
            return NULL;
        return boost::python::incref(loglevel_map[self->impl_.level()].ptr());
    }

    static PyObject* emit(PyObject* _self, PyObject* _record)
    {
        CppLoggerHandler* const self(get_self(_self));
        if (!self)
            return NULL;
        boost::python::object record(
            boost::python::borrowed(_record));

        try
        {
            boost::python::handle<> _level(
                    PyObject_GetAttrString(_record, "levelno"));
            boost::python::object msg(
                boost::python::getattr(
                    boost::python::object(
                        boost::python::borrowed(_self)), "format")(record));
            self->impl_.log(translate_level(_level.get()), "%s",
                    boost::python::extract<char const*>(msg)());
        }
        catch (boost::python::error_already_set const&)
        {
            return NULL;
        }
        return boost::python::incref(Py_None);
    }

    static PyObject* flush(PyObject* _self, PyObject* arg)
    {
        CppLoggerHandler* const self(get_self(_self));
        if (!self)
            return NULL;
        self->impl_.flush();
        return boost::python::incref(Py_None);
    }

    static PyObject* close(PyObject* _self, PyObject* arg)
    {
        CppLoggerHandler* const self(get_self(_self));
        if (!self)
            return NULL;
        return boost::python::incref(Py_None);
    }

    static PyObject* translateLevelValue(PyObject* _self, PyObject* arg)
    try
    {
        return boost::python::incref(boost::python::object(translate_level(arg)).ptr());
    }
    catch (boost::python::error_already_set const&)
    {
        return NULL;
    }

    static PyObject* get_logger(PyObject* _self, void* context)
    {
        CppLoggerHandler* const self(get_self(_self));
        if (!self)
            return NULL;
        return boost::python::reference_existing_object::apply<Logger*>::type()(self->impl_);
    }


#ifndef HAVE_OLD_FASHIONED_LOGGER_CLASSES
    static PyObject* get___dict__(PyObject *_self, void* context)
    {
        CppLoggerHandler* const self(get_self(_self));
        if (!self)
            return NULL;
        return boost::python::incref(self->__dict__.get());
    }
#endif

    CppLoggerHandler(Logger& impl)
        :
#ifndef HAVE_OLD_FASHIONED_LOGGER_CLASSES
        __weakreflist__(0),
        __dict__(PyDict_New()),
#endif
        impl_(impl)
    {}

protected:
    static enum Logger::level translate_level(PyObject* _level)
    {
        enum Logger::level retval(Logger::L_OFF);
        boost::python::object level(boost::python::borrowed(_level));
        boost::python::object closest;

        BOOST_FOREACH (loglevel_map_type::value_type const& i, loglevel_map)
        {
            if (i.second <= level && closest < i.second)
            {
                retval = i.first;
                closest = i.second;
            }
        }
        return retval;
    }


protected:
#ifdef HAVE_OLD_FASHIONED_LOGGER_CLASSES
    static boost::python::handle<> __class__;
#else
    PyObject_VAR_HEAD
    static PyTypeObject __class__;
    static std::string __name__;
    PyObject *__weakreflist__;
    boost::python::handle<> __dict__;
#endif
    static PyMethodDef __methods__[];
    static PyGetSetDef __getsets__[];
    Logger& impl_;
};

#ifndef HAVE_OLD_FASHIONED_LOGGER_CLASSES
std::string CppLoggerHandler::__name__;
#endif

PyMethodDef CppLoggerHandler::__methods__[] = {
#ifdef HAVE_OLD_FASHIONED_LOGGER_CLASSES
    { "__init__", (PyCFunction)CppLoggerHandler::__init__, METH_O, "" },
#endif
    { "createLock", (PyCFunction)CppLoggerHandler::createLock, METH_NOARGS, "" },
    { "acquire", (PyCFunction)CppLoggerHandler::acquire, METH_NOARGS, "" },
    { "release", (PyCFunction)CppLoggerHandler::release, METH_NOARGS, "" },
    { "setLevel", (PyCFunction)CppLoggerHandler::setLevel, METH_O, "" },
    { "getLevel", (PyCFunction)CppLoggerHandler::getLevel, METH_NOARGS, "" },
    { "emit", (PyCFunction)CppLoggerHandler::emit, METH_O, "" },
    { "flush", (PyCFunction)CppLoggerHandler::flush, METH_NOARGS, "" },
    { "close", (PyCFunction)CppLoggerHandler::close, METH_NOARGS, "" },
    { "translateLevelValue", (PyCFunction)CppLoggerHandler::translateLevelValue, METH_O | METH_STATIC, "" },
    { NULL, NULL }
};

PyGetSetDef CppLoggerHandler::__getsets__[] = {
#ifndef HAVE_OLD_FASHIONED_LOGGER_CLASSES
    { const_cast<char*>("__dict__"), (getter)&CppLoggerHandler::get___dict__, NULL, NULL },
#endif
    { const_cast<char*>("logger"), (getter)&CppLoggerHandler::get_logger, NULL, NULL },
    { NULL, NULL }
};

#ifdef HAVE_OLD_FASHIONED_LOGGER_CLASSES
boost::python::handle<> CppLoggerHandler::__class__;
#else
PyTypeObject CppLoggerHandler::__class__ = {
	PyObject_HEAD_INIT(&PyType_Type)
	0,					/* ob_size */
	0,                  /* tp_name */
	sizeof(CppLoggerHandler), /* tp_basicsize */
	0,					/* tp_itemsize */
	/* methods */
	(destructor)&CppLoggerHandler::__dealloc__, /* tp_dealloc */
	0,					/* tp_print */
	0,					/* tp_getattr */
	0,					/* tp_setattr */
	0,					/* tp_compare */
	0,					/* tp_repr */
	0,					/* tp_as_number */
	0,                  /* tp_as_sequence */
	0,					/* tp_as_mapping */
	0,					/* tp_hash */
	0,					/* tp_call */
	0,					/* tp_str */
	PyObject_GenericGetAttr,		/* tp_getattro */
	0,					/* tp_setattro */
	0,					/* tp_as_buffer */
	Py_TPFLAGS_HAVE_CLASS | Py_TPFLAGS_HAVE_WEAKREFS | Py_TPFLAGS_BASETYPE, /* tp_flags */
	0,					/* tp_doc */
	CppLoggerHandler::__traverse__,              	/* tp_traverse */
	0,					/* tp_clear */
	0,                  /* tp_richcompare */
	offsetof(CppLoggerHandler, __weakreflist__),    /* tp_weaklistoffset */
	0,                  /* tp_iter */
	0,                  /* tp_iternext */
	CppLoggerHandler::__methods__,		        	/* tp_methods */
	0,					/* tp_members */
    CppLoggerHandler::__getsets__, /* tp_getset */
    0,                  /* tp_base */
    0,                  /* tp_dict */
    0,                  /* tp_descr_get */
    0,                  /* tp_descr_set */
    offsetof(CppLoggerHandler, __dict__),           /* tp_dictoffset */
    CppLoggerHandler::__init__, /* tp_init */
    0,                  /* tp_alloc */
    CppLoggerHandler::__new__,  /*tp_new */
    0                   /* tp_free */
};
#endif

boost::python::object
register_logger_handler_class(char const* name)
{
    using namespace boost::python;

    import_logging_module();

    PyObject* klass(
        CppLoggerHandler::__class_init__(
            name, reinterpret_cast<PyObject*>(scope().ptr())));
    boost::python::object retval(borrowed(klass));
    scope().attr(name) = retval;
    return retval;
}

boost::python::objects::class_base
register_python_appender_class(char const* name)
{
    using namespace boost::python;
    typedef PythonAppender impl_type;

    import_logging_module();

    return class_<impl_type, bases<LogAppender>,
                  boost::shared_ptr<impl_type>, boost::noncopyable>(
                  name, init<boost::python::object>())
        ;
}

} // namespace binding
