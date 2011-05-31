#ifndef BINDING_REACTION_RECORDER_CONVERTER_HPP
#define BINDING_REACTION_RECORDER_CONVERTER_HPP

#include <boost/python.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <boost/scoped_ptr.hpp>
#include "peer/utils.hpp"

namespace binding {

template<typename Tbase_>
class ReactionRecorderWrapper: public Tbase_
{
public:
    typedef Tbase_ wrapped_type;
    typedef typename wrapped_type::reaction_record_type reaction_record_type;

public:
    virtual ~ReactionRecorderWrapper() {}

    virtual void operator()(reaction_record_type const& rr)
    {
        boost::python::decref(PyObject_CallObject(callable_, boost::python::make_tuple(boost::python::object(rr)).ptr()));
    }

    ReactionRecorderWrapper(PyObject* callable): callable_(callable) {}

private:
    PyObject* callable_;
};

template<typename Tbase_>
struct reaction_recorder_lvalue_converter
{
    typedef ReactionRecorderWrapper<Tbase_> native_type;

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

template<typename Tbase_>
struct reaction_recorder_converter
{
    typedef boost::shared_ptr<Tbase_> native_type;

    static void* convertible(PyObject* ptr)
    {
        if (!PyCallable_Check(ptr))
        {
            return NULL;
        }

        return ptr;
    }
    
    static void construct(PyObject* ptr,
                          boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        data->stage1.convertible =
            new(data->storage.bytes) native_type(
                new ReactionRecorderWrapper<Tbase_>(
                    static_cast<PyObject*>(data->stage1.convertible)));
    }
};

template<typename Timpl_>
void register_reaction_recorder_converter()
{
    using namespace boost::python;
    typedef Timpl_ impl_type;

    peer::util::to_native_lvalue_converter<impl_type, reaction_recorder_lvalue_converter<impl_type> >();
    peer::util::to_native_converter<boost::shared_ptr<impl_type>, reaction_recorder_converter<impl_type> >();
}

} // namespace binding

#endif /* BINDING_REACTION_RECORDER_CONVERTER_HPP */
