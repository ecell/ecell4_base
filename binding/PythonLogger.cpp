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

    virtual Logger* operator()(char const* name) const;

    virtual char const* get_name() const;
};

// XXX: logging.Handler is a old-style class... *sigh*
class LoggingHandler: public PyInstanceObject
{
public:
    void* operator new(size_t)
    {
        PyInstanceObject *inst;

        boost::python::handle<> dict(PyDict_New());
        inst = reinterpret_cast<PyInstanceObject*>(
            _PyObject_GC_Malloc(sizeof(LoggingHandler)));
        if (inst)
        {
            inst = reinterpret_cast<PyInstanceObject*>(
                PyObject_INIT(reinterpret_cast<PyObject*>(inst),
                                 &PyInstance_Type));
        }
        if (!inst)
            throw std::bad_alloc();
        inst->in_weakreflist = NULL;
        inst->in_class = reinterpret_cast<PyClassObject*>(
            boost::python::incref(__real_class__.get()));
        inst->in_dict = boost::python::incref(dict.get());
        _PyObject_GC_TRACK(inst);
        return inst;
    }

    void operator delete(void* ptr)
    {
        PyInstance_Type.tp_dealloc(reinterpret_cast<PyObject*>(ptr));
    }

    static PyObject* __class_init__(const char* name, PyObject* mod)
    {
        using namespace boost::python;

        import_logging_module();

        __name__ = static_cast<std::string>(
            extract<std::string>(object(borrowed(mod)).attr("__name__")))
            + "." + name;

        __type__ = PyClass_Type;
        __type__.tp_call = (ternaryfunc)&__new__;

        boost::python::handle<> dict(PyDict_New());
        if (PyDict_SetItemString(dict.get(), "__module__", mod))
            return NULL;

        boost::python::object bases(
            boost::python::make_tuple(
                getattr(logging_module, "Handler")));
        boost::python::object _name(__name__);

        boost::python::handle<> klass(
            PyClass_New(bases.ptr(), dict.get(), _name.ptr()));
        if (PyErr_Occurred())
            return NULL;

        boost::python::handle<> real_klass(
            PyClass_New(bases.ptr(), dict.get(), _name.ptr()));
        if (PyErr_Occurred())
            return NULL;

        klass->ob_type = &__type__;

        for (PyMethodDef* meth = __methods__; meth->ml_name; ++meth)
        {
            boost::python::handle<> _meth(PyMethod_New(
                PyCFunction_New(meth, real_klass.get()), NULL, real_klass.get()));
            if (PyDict_SetItemString(dict.get(), meth->ml_name, _meth.get()))
            {
                return NULL;
            }
        }

        __class__ = klass;
        __real_class__ = real_klass;

        return klass.get();
    }

    static PyObject* __new__(PyClassObject* klass, PyObject* arg, PyObject* kwarg)
    {
        if (PyTuple_Size(arg) != 1)
        {
            PyErr_SetString(PyExc_TypeError, "the number of arguments must be 1");
            return NULL;
        }
        boost::python::object name(
                boost::python::borrowed(PyTuple_GetItem(arg, 0)));

        PyObject* const retval(reinterpret_cast<PyObject*>(
            new LoggingHandler(
                Logger::get_logger(
                    boost::python::extract<char const*>(name)))));

        boost::python::handle<> bases(
            boost::python::borrowed(klass->cl_bases));
   
        if (!PyTuple_Check(bases.get()))
            return NULL;

        boost::python::handle<> super_arg(PyTuple_Pack(1, retval));

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
                Py_XDECREF(retval);
                return NULL; 
            }
        }

        return retval;
    }

    static void __dealloc__(LoggingHandler* self)
    {
        delete self;
    }

    static PyObject* createLock(LoggingHandler* self)
    {
        return boost::python::incref(Py_None);
    }

    static PyObject* acquire(LoggingHandler* self)
    {
        return boost::python::incref(Py_None);
    }

    static PyObject* release(LoggingHandler* self)
    {
        return boost::python::incref(Py_None);
    }

    static PyObject* setLevel(LoggingHandler* self, PyObject* level)
    {
        boost::python::object _level(boost::python::borrowed(level));
        boost::python::object closest;
        BOOST_FOREACH (loglevel_map_type::value_type const& i, loglevel_map)
        {
            if (i.second <= _level && closest < i.second)
            {
                closest = i.second;
            }
        }
        if (closest.ptr() == Py_None)
        {
            PyErr_SetString(PyExc_ValueError, "invalid loglevel");
            return NULL;
        }
        return boost::python::incref(Py_None);
    }

    static PyObject* emit(LoggingHandler* self, PyObject* record)
    {
        boost::python::object msg(
            boost::python::getattr(
                boost::python::object(
                    boost::python::borrowed(
                        reinterpret_cast<PyObject*>(self))), "format")(
                            boost::python::object(
                                boost::python::borrowed(record))));
        self->impl_.log(self->level_, "%s",
                boost::python::extract<char const*>(msg)());
        return boost::python::incref(Py_None);
    }

    static PyObject* flush(LoggingHandler* self)
    {
        return boost::python::incref(Py_None);
    }

    static PyObject* close(LoggingHandler* self)
    {
        return boost::python::incref(Py_None);
    }

    LoggingHandler(Logger& impl): impl_(impl), level_(Logger::L_INFO) {}

protected:
    static PyTypeObject __type__;
    static boost::python::handle<> __class__;
    static boost::python::handle<> __real_class__;
    static PyMethodDef __methods__[];
    static std::string __name__;
    Logger& impl_;
    enum Logger::level level_;
};

std::string LoggingHandler::__name__;

PyTypeObject LoggingHandler::__type__;

boost::python::handle<> LoggingHandler::__class__;
boost::python::handle<> LoggingHandler::__real_class__;

PyMethodDef LoggingHandler::__methods__[] = {
    { "createLock", (PyCFunction)LoggingHandler::createLock, METH_NOARGS, "" },
    { "acquire", (PyCFunction)LoggingHandler::acquire, METH_NOARGS, "" },
    { "release", (PyCFunction)LoggingHandler::release, METH_NOARGS, "" },
    { "setLevel", (PyCFunction)LoggingHandler::setLevel, METH_O, "" },
    { "emit", (PyCFunction)LoggingHandler::emit, METH_O, "" },
    { "flush", (PyCFunction)LoggingHandler::flush, METH_O, "" },
    { "close", (PyCFunction)LoggingHandler::close, METH_O, "" },
    { NULL, NULL }
};

PythonLoggerFactory::~PythonLoggerFactory()
{
}

Logger* PythonLoggerFactory::operator()(char const* name) const
{
    using namespace boost::python;
    object getLogger(getattr(logging_module, "getLogger"));
    return new PythonLogger(getLogger(name));
}

char const* PythonLoggerFactory::get_name() const
{
    return "PythonLogger";
}

boost::python::objects::class_base
register_logger_factory_class(char const* name)
{
    using namespace boost::python;
    typedef LoggerFactory impl_type;

    return class_<impl_type, boost::shared_ptr<impl_type>, boost::noncopyable>(name, no_init)
        .add_property("name", &impl_type::get_name)
        .def("register_logger_factory", &impl_type::register_logger_factory)
        .staticmethod("register_logger_factory")
        ;
}

boost::python::object
register_logger_handler_class(char const* name)
{
    using namespace boost::python;
    PyObject* klass(
        LoggingHandler::__class_init__(
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
