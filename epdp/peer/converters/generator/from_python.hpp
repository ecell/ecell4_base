#ifndef PEER_CONVERTERS_GENERATOR_TO_PYTHON_HPP
#define PEER_CONVERTERS_GENERATOR_TO_PYTHON_HPP

namespace peer { namespace converters {

struct pyiterator_to_generator_converter
{
    typedef Twrapper_ native_type;
    typedef Tgen_ generator_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PyObject_TypeCheck(pyo, &native_type::__class__))
        {
            return 0;
        }
        return pyo;
    }

    static void construct(PyObject* pyo, 
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        void* storage(reinterpret_cast<
            boost::python::converter::rvalue_from_python_storage<Tgen_>* >(
                data)->storage.bytes);
        new (storage) Tgen_(reinterpret_cast<native_type*>(pyo)->impl_);
        data->convertible = storage;
    }
};

} // namespace peer::converters

#endif /* PEER_CONVERTERS_GENERATOR_TO_PYTHON_HPP */
