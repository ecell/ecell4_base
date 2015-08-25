from libcpp.string cimport string

cdef extern from "ecell4/core/extras.hpp" namespace "ecell4::extras":
    string load_version_information(string&)
