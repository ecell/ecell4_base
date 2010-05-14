#ifndef PEER_CONVERTERS_GENERATOR_TO_PYTHON_HPP
#define PEER_CONVERTERS_GENERATOR_TO_PYTHON_HPP

namespace peer { namespace converters {

template<typename Tgen_>
struct generator_to_pyiterator_converter 
{
    typedef peer::wrappers::generator_wrapper<Tgen_> native_type;
    typedef Tgen_ generator_type;

    static PyObject* convert(Tgen_ const& impl)
    {
        return native_type::create(impl);
    }

    static PyTypeObject* get_pytype()
    {
        return &native_type::__class__;
    }
};

} } // namespace peer::converters

#endif /* PEER_CONVERTERS_GENERATOR_TO_PYTHON_HPP */
