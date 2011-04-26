#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/python.hpp>
#include <boost/date_time/posix_time/conversion.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include "peer/util/to_native_converter.hpp"
#include "binding_common.hpp"

namespace binding {

static boost::python::object datetime_module;

static void import_datetime_module()
{
    if (datetime_module.ptr() != Py_None)
        return;

    datetime_module = boost::python::object(boost::python::borrowed(PyImport_Import(boost::python::handle<>(PyString_FromString("datetime")).get())));
    if (PyErr_Occurred())
    {
        boost::python::throw_error_already_set();
    }
}

template<typename Timpl>
static void log_appender_call(Timpl& self, enum Logger::level lv,
                              char const* name, PyObject* _chunks)
{
    using namespace boost::python;

    handle<> iter(allow_null(PyObject_GetIter(_chunks)));
    if (!iter)
    {
        PyErr_SetString(PyExc_TypeError, "chunks must be iterable");
        throw_error_already_set();    
    }

    std::vector<char const*> chunks;
    for (;;)
    {
        handle<> item(allow_null(PyIter_Next(iter.get())));
        if (!item)
            break;
        chunks.push_back(extract<char const*>(item.get()));
    }

    chunks.push_back(0);
    self(lv, name, &chunks.front());
}

struct posix_timer_converter
{
    typedef boost::posix_time::ptime native_type;

    static void* convertible(PyObject* ptr)
    {
        return ptr;
    }
    
    static void construct(PyObject* ptr,
                          boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        if (PyInt_Check(ptr))
        {
            data->stage1.convertible = new(data->storage.bytes) native_type(
                boost::posix_time::from_time_t(PyInt_AsLong(ptr)));
        }
        else if (PyObject_IsSubclass(reinterpret_cast<PyObject*>(Py_TYPE(ptr)),
                boost::python::getattr(datetime_module, "datetime").ptr()))
        {
            boost::python::object _ptr(boost::python::borrowed(ptr));
            data->stage1.convertible = new(data->storage.bytes) native_type(
                boost::gregorian::date(
                    boost::python::extract<int>(
                        boost::python::getattr(_ptr, "year"))(),
                    boost::python::extract<int>(
                        boost::python::getattr(_ptr, "month"))(),
                    boost::python::extract<int>(
                        boost::python::getattr(_ptr, "day"))()),
                boost::posix_time::time_duration(
                    boost::python::extract<int>(
                        boost::python::getattr(_ptr, "hour"))(),
                    boost::python::extract<int>(
                        boost::python::getattr(_ptr, "minute"))(),
                    boost::python::extract<int>(
                        boost::python::getattr(_ptr, "second"))(),
                    boost::python::extract<int>(
                        boost::python::getattr(_ptr, "microsecond"))() * 1000));
        }
        else
        {
            PyErr_SetString(PyExc_TypeError, "argument must be either in or datetime");
            boost::python::throw_error_already_set();
        }
    }
};

boost::python::objects::class_base
register_log_appender_class(char const* name)
{
    using namespace boost::python;
    typedef LogAppender impl_type;

    import_datetime_module();

    peer::util::to_native_converter<boost::posix_time::ptime, posix_timer_converter>();

    return class_<impl_type, boost::shared_ptr<impl_type>, boost::noncopyable>(name, no_init)
        .def("flush", &impl_type::flush)
        .def("__call__", &log_appender_call<impl_type>)
        ;
}

} // namespace binding
