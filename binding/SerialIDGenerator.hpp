#ifndef PEER_SERIAL_ID_GENERATOR_HPP
#define PEER_SERIAL_ID_GENERATOR_HPP

#include <boost/python.hpp>
#include "../SerialIDGenerator.hpp"

template<typename T>
inline void register_serial_id_generator_class(char const* class_name)
{
    using namespace boost::python;

    class_<SerialIDGenerator<T> >(class_name, init<int>())
        .def("__call__", &SerialIDGenerator<T>::operator())
        ;
}

#endif /* PEER_SERIAL_ID_GENERATOR_HPP */
