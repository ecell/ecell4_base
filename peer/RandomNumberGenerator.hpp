#ifndef PEER_RANDOMNUMBERGENERATOR_HPP
#define PEER_RANDOMNUMBERGENERATOR_HPP

#include <boost/python.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include "utils/reference_or_instance.hpp"

template<typename Trng_>
struct RandomNumberGenerator
{
    PyObject* normal(Real loc, Real scale, boost::python::object size)
    {
        boost::normal_distribution<Real> dist(loc, scale);
        if (size.ptr() == Py_None) {
            return boost::python::incref(boost::python::object(dist(*this)).ptr());
        } else {
            ssize_t len = 0;
            std::size_t num_samples = 1;
            PyObject* retval = 0;
            len = PySequence_Size(size.ptr());
            // cancel the error that would occur if size wasn't be a sequence
            if (PyErr_Occurred())
            {
                PyErr_Clear();
            }
            if (len == -1 && PyNumber_Check(size.ptr()))
            {
                npy_intp dims[1] = { PyNumber_AsSsize_t(size.ptr(), NULL) };
                if (dims[0] < 0)
                {
                    PyErr_SetString(PyExc_ValueError, "expected size >= 0");
                    boost::python::throw_error_already_set();
                }
                num_samples = dims[0];
                retval = PyArray_New(
                        &PyArray_Type, 1, dims,
                        peer::util::get_numpy_typecode<Real>::value,
                        NULL, NULL, 0, NPY_CARRAY, NULL);
            }
            else if (len >= 0)
            {
                if (len > NPY_MAXDIMS)
                {
                    PyErr_Format(PyExc_ValueError, "sequence too large; must be smaller than %d", NPY_MAXDIMS);
                    boost::python::throw_error_already_set();
                }
                npy_intp dims[len];
                for (ssize_t i = 0; i < len; ++i)
                {
                    dims[i] = boost::python::extract<npy_intp>(size[i]);
                    num_samples *= dims[i];
                }
                retval = PyArray_New(
                        &PyArray_Type, len, dims,
                        peer::util::get_numpy_typecode<Real>::value,
                        NULL, NULL, 0, NPY_CARRAY, NULL);
            }
            Real* data(reinterpret_cast<Real*>(PyArray_DATA(retval)));
            for (std::size_t i = 0; i < num_samples; ++i)
            {
                data[i] = dist(*this);
            }
            return retval;
        }
    }

    Real uniform(Real min, Real max)
    {
        return boost::uniform_real<Real>(min, max)(static_cast<Trng_&>(rng_));
    }

    int uniform_int(int min, int max)
    {
        return boost::uniform_int<int>(min, max)(static_cast<Trng_&>(rng_));
    }

    Real operator()()
    {
        return uniform(0.0, 1.0);
    }

    void seed(typename Trng_::result_type val)
    {
        static_cast<Trng_&>(rng_).seed(val);
    }

    RandomNumberGenerator(Trng_& rng): rng_(rng) {}

    RandomNumberGenerator(): rng_() {}

    static void __register_class(const char* name)
    {
        using namespace boost::python;
        class_<RandomNumberGenerator, boost::noncopyable>(name, init<>())
            .def("normal", &RandomNumberGenerator::normal)
            .def("seed", &RandomNumberGenerator::seed)
            .def("uniform", &RandomNumberGenerator::uniform)
            .def("uniform_int", &RandomNumberGenerator::uniform_int)
            .def("__call__", &RandomNumberGenerator::operator())
            ;
    }

    reference_or_instance<Trng_> rng_;
};

#endif /* PEER_RANDOMNUMBERGENERATOR_HPP */
