#ifndef PYLEX_HELPER_HPP
#define PYLEX_HELPER_HPP

#include <Python.h>
#include <string>
#include <iostream>
#include <sstream>

namespace pyrex_helper {

template<typename T_>
inline void ref_set(T_& ref, const T_& val) { ref = val; }

template<typename T_>
inline PyObject* pystr_from_repr(const T_* obj)
{
    std::ostringstream s;
    s << (*obj);
    std::string res(s.str());
    return PyString_FromStringAndSize(res.data(), res.size());
}

} // namespace pyrex_helper

#endif /* PYLEX_HELPER_HPP */
