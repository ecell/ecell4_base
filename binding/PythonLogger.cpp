#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <map>
#include <cstdio>
#include <boost/foreach.hpp>
#include "PythonLogger.hpp"
#include "binding_common.hpp"

namespace binding {

typedef std::map<enum Logger::level, boost::python::object> loglevel_map_type;

static boost::python::object logging_module;
static loglevel_map_type loglevel_map;

static void import_logging_module()
{
    if (logging_module.ptr() != Py_None)
        return;

    logging_module = boost::python::object(boost::python::borrowed(PyImport_Import(PyString_FromString("logging"))));
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

class PythonLogger: public Logger
{
public:
    virtual ~PythonLogger() {}

    virtual void level(enum level level)
    {
        logger_.attr("setLevel")(loglevel_map[level]);
    }

    virtual enum level level() const
    {
        return Logger::L_DEBUG;
    }

    virtual void logv(enum level lv, char const* format, va_list ap)
    {
        char buf[2048];
        std::vsnprintf(buf, sizeof(buf), format, ap);
        logger_.attr("log")(loglevel_map[lv], "%s", buf);
    }

    virtual void flush()
    {
        // do nothing for now.
    }

    PythonLogger(boost::python::object logger): logger_(logger) {}

private:
    boost::python::object logger_;
};

class PythonLoggerFactory: public LoggerFactory
{
public:
    virtual ~PythonLoggerFactory();

    virtual void level(enum Logger::level);

    virtual Logger* operator()(char const* name) const;

    virtual char const* get_name() const;

    PythonLoggerFactory();

private:
    enum Logger::level level_;
    boost::python::object getLogger;
};

// XXX: logging.Handler is a old-style class... *sigh*
class CppLoggerHandler
{
public:
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
                    PyDescr_NewMethod(&PyInstance_Type, meth)));
            if (!_meth)
                return NULL;
            if (PyDict_SetItemString(dict.get(), meth->ml_name, _meth.get()))
                return NULL;
        }

        __class__ = klass;

        return klass.get();
    }

    static void __dealloc__(void* ptr)
    {
        delete reinterpret_cast<CppLoggerHandler*>(ptr);
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

    static PyObject* __init__(PyObject* _self, PyObject* name)
    {
        BOOST_ASSERT(PyInstance_Check(_self));

        CppLoggerHandler* self(new CppLoggerHandler(
            Logger::get_logger(
                boost::python::extract<char const*>(name))));
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
        enum Logger::level const level(get_level(_level));
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
            self->impl_.log(get_level(_level.get()), "%s",
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

    CppLoggerHandler(Logger& impl): impl_(impl) {}

protected:
    static enum Logger::level get_level(PyObject* _level)
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
    static boost::python::handle<> __class__;
    static PyMethodDef __methods__[];
    Logger& impl_;
};

boost::python::handle<> CppLoggerHandler::__class__;

PyMethodDef CppLoggerHandler::__methods__[] = {
    { "__init__", (PyCFunction)CppLoggerHandler::__init__, METH_O, "" },
    { "createLock", (PyCFunction)CppLoggerHandler::createLock, METH_NOARGS, "" },
    { "acquire", (PyCFunction)CppLoggerHandler::acquire, METH_NOARGS, "" },
    { "release", (PyCFunction)CppLoggerHandler::release, METH_NOARGS, "" },
    { "setLevel", (PyCFunction)CppLoggerHandler::setLevel, METH_O, "" },
    { "getLevel", (PyCFunction)CppLoggerHandler::getLevel, METH_NOARGS, "" },
    { "emit", (PyCFunction)CppLoggerHandler::emit, METH_O, "" },
    { "flush", (PyCFunction)CppLoggerHandler::flush, METH_NOARGS, "" },
    { "close", (PyCFunction)CppLoggerHandler::close, METH_NOARGS, "" },
    { NULL, NULL }
};

PythonLoggerFactory::~PythonLoggerFactory()
{
}

void PythonLoggerFactory::level(enum Logger::level level)
{
    level_ = level;
}

Logger* PythonLoggerFactory::operator()(char const* name) const
{
    using namespace boost::python;
    object logger(getLogger(name));
    logger.attr("setLevel")(loglevel_map[level_]);
    return new PythonLogger(logger);
}

char const* PythonLoggerFactory::get_name() const
{
    return "PythonLogger";
}

PythonLoggerFactory::PythonLoggerFactory()
    : level_(Logger::L_DEBUG), getLogger(getattr(logging_module, "getLogger"))
{
}

static PyObject* dummy_getter(PyObject*)
{
    PyErr_SetString(PyExc_AttributeError, "try to read a write-only property");
    return NULL;
}

boost::python::objects::class_base
register_logger_factory_class(char const* name)
{
    using namespace boost::python;
    typedef LoggerFactory impl_type;

    return class_<impl_type, boost::shared_ptr<impl_type>, boost::noncopyable>(name, no_init)
        .add_property("level", &dummy_getter, &impl_type::level)
        .add_property("name", &impl_type::get_name)
        .def("register_logger_factory", &impl_type::register_logger_factory)
        .staticmethod("register_logger_factory")
        .def("get_logger_factory", &impl_type::get_logger_factory)
        .staticmethod("get_logger_factory")
        ;
}

boost::python::objects::enum_base
register_logger_level_enum(char const* name)
{
    using namespace boost::python;
    return enum_<enum   Logger::level>(name)
        .value("OFF",     Logger::L_OFF)
        .value("DEBUG",   Logger::L_DEBUG)
        .value("INFO",    Logger::L_INFO)
        .value("WARNING", Logger::L_WARNING)
        .value("ERROR",   Logger::L_ERROR)
        .value("FATAL",   Logger::L_FATAL)
        ;
}

boost::python::object
register_logger_handler_class(char const* name)
{
    using namespace boost::python;
    PyObject* klass(
        CppLoggerHandler::__class_init__(
            name, reinterpret_cast<PyObject*>(scope().ptr())));
    boost::python::object retval(borrowed(klass));
    scope().attr(name) = retval;
    return retval;
}

boost::python::objects::class_base
register_null_logger_factory_class(char const* name)
{
    using namespace boost::python;
    typedef NullLoggerFactory impl_type;

    return class_<impl_type, bases<LoggerFactory>, boost::shared_ptr<impl_type>, boost::noncopyable>(name)
        ;
}

boost::python::objects::class_base
register_python_logger_factory_class(char const* name)
{
    using namespace boost::python;
    typedef PythonLoggerFactory impl_type;

    import_logging_module();

    return class_<impl_type, bases<LoggerFactory>, boost::shared_ptr<impl_type>, boost::noncopyable>(name)
        ;
}

} // namespace binding
