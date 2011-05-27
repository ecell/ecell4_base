#ifndef BINDING_VOLUME_CLEARER_CONVERTER_HPP
#define BINDING_VOLUME_CLEARER_CONVERTER_HPP

#include <boost/python.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <boost/scoped_ptr.hpp>
#include "peer/utils.hpp"

namespace binding {

template<typename Tbase_>
class VolumeClearerWrapper: public Tbase_
{
public:
    typedef Tbase_ wrapped_type;
    typedef typename wrapped_type::particle_shape_type particle_shape_type;
    typedef typename wrapped_type::particle_id_type particle_id_type;

public:
    virtual ~VolumeClearerWrapper() {}

    virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore)
    {
        PyObject* retobj(PyObject_CallObject(callable_, boost::python::make_tuple(boost::python::object(shape), boost::python::object(ignore)).ptr()));
        bool const retval(retobj == Py_True);
        boost::python::decref(retobj);
        return retval;
    }

    virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore0, particle_id_type const& ignore1)
    {
        PyObject* retobj(PyObject_CallObject(callable_, boost::python::make_tuple(boost::python::object(shape), boost::python::object(ignore0), boost::python::object(ignore1)).ptr()));
        bool const retval(retobj == Py_True);
        boost::python::decref(retobj);
        return retval;
    }

    VolumeClearerWrapper(PyObject* callable): callable_(callable) {}

private:
    PyObject* callable_;
};

template<typename Tbase_>
struct volume_clearer_converter
{
    typedef VolumeClearerWrapper<Tbase_> native_type;

    static PyTypeObject const* expected_pytype()
    {
        return &PyBaseObject_Type;
    }

    static void* convert(PyObject* pyo)
    {
        if (!pyo || pyo == Py_None)
        {
            return 0;
        }

        if (!PyCallable_Check(pyo))
        {
            boost::python::throw_error_already_set();
        }

        native_type* retval(new native_type(pyo));
        peer::util::install_instance_holder<boost::scoped_ptr<native_type> >(pyo, boost::in_place(retval));
        return retval;
    }
};

template<typename Timpl_>
void register_volume_clearer_converter()
{
    using namespace boost::python;
    typedef Timpl_ impl_type;

    peer::util::to_native_lvalue_converter<impl_type, volume_clearer_converter<impl_type> >();
}

} // namespace binding

#endif /* BINDING_VOLUME_CLEARER_CONVERTER_HPP */
