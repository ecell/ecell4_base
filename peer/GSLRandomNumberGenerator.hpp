#ifndef PEER_GSLRANDOMNUMBERGENERATOR_HPP
#define PEER_GSLRANDOMNUMBERGENERATOR_HPP

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace peer {

struct GSLRandomNumberGenerator
{
    typedef boost::shared_ptr<gsl_rng> rng_handle;

    virtual ~GSLRandomNumberGenerator() {}

    PyObject* normal(Real loc, Real scale, boost::python::object size)
    {
        if (size.ptr() == Py_None) {
            return boost::python::incref(boost::python::object(
                gsl_ran_gaussian(rng_.get(), scale) + loc).ptr());
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
                data[i] = gsl_ran_gaussian(rng_.get(), scale) + loc;
            }
            return retval;
        }
    }

    Real uniform(Real min, Real max)
    {
        return gsl_rng_uniform(rng_.get()) * (max - min) + min;
    }

    int uniform_int(int min, int max)
    {
        return gsl_rng_uniform_int(rng_.get(), max - min + 1) + min;
    }

    Real operator()()
    {
        return gsl_rng_uniform(rng_.get());
    }

    void seed(unsigned long int val)
    {
        gsl_rng_set(rng_.get(), val);
    }

    GSLRandomNumberGenerator(gsl_rng* rng): rng_(rng, gsl_rng_free) {}

    template<gsl_rng_type const*& Prng_>
    static GSLRandomNumberGenerator create()
    {
        return GSLRandomNumberGenerator(gsl_rng_alloc(Prng_));
    }

    template<gsl_rng_type const*& Prng_>
    static void __register_class(const char* name)
    {
        using namespace boost::python;
        class_<GSLRandomNumberGenerator>(name, no_init)
            .def("create", &GSLRandomNumberGenerator::create<Prng_>,
                           return_value_policy<return_by_value>())
            .staticmethod("create")
            .def("normal", &GSLRandomNumberGenerator::normal)
            .def("seed", &GSLRandomNumberGenerator::seed)
            .def("uniform", &GSLRandomNumberGenerator::uniform)
            .def("uniform_int", &GSLRandomNumberGenerator::uniform_int)
            .def("__call__", &GSLRandomNumberGenerator::operator())
            ;

    }

    rng_handle rng_;
};

} // namespace peer

#endif /* PEER_GSLRANDOMNUMBERGENERATOR_HPP */
