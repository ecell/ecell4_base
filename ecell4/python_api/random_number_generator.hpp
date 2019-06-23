#ifndef ECELL4_PYTHON_API_RANDOM_NUMBER_GENERATOR_HPP
#define ECELL4_PYTHON_API_RANDOM_NUMBER_GENERATOR_HPP

#include <pybind11/pybind11.h>

namespace ecell4
{

namespace python_api
{

    template<class Base = ecell4::RandomNumberGenerator>
    class PyRandomNumberGenerator: public Base
    {
    public:
        using Base::Base;

        Real random()
        {
            PYBIND11_OVERLOAD_PURE(Real, Base, random,);
        }

        Real uniform(Real min, Real max)
        {
            PYBIND11_OVERLOAD_PURE(Real, Base, uniform, min, max);
        }

        Integer uniform_int(Integer min, Integer max)
        {
            PYBIND11_OVERLOAD_PURE(Integer, Base, uniform_int, min, max);
        }

        Real gaussian(Real sigma, Real mean = 0.0)
        {
            PYBIND11_OVERLOAD_PURE(Real, Base, gaussian, sigma, mean);
        }

        Integer binomial(Real p, Integer n)
        {
            PYBIND11_OVERLOAD_PURE(Integer, Base, binomial, p, n);
        }

        Real3 direction3d(Real length = 1.0)
        {
            PYBIND11_OVERLOAD_PURE(Real3, Base, direction3d, length);
        }

        void seed(Integer val)
        {
            PYBIND11_OVERLOAD_PURE(void, Base, seed, val);
        }

        void seed()
        {
            PYBIND11_OVERLOAD_PURE(void, Base, seed,);
        }

#ifdef WITH_HDF5
        void save(H5::H5Location* root) const
        {
            PYBIND11_OVERLOAD_PURE(void, Base, save, root);
        }

        void load(const H5::H5Location& root)
        {
            PYBIND11_OVERLOAD_PURE(void, Base, load, root);
        }

        void save(const std::string& filename) const
        {
            PYBIND11_OVERLOAD_PURE(void, Base, save, filename);
        }

        void load(const std::string& filename)
        {
            PYBIND11_OVERLOAD_PURE(void, Base, load, filename);
        }
#else
        void save(const std::string& filename) const
        {
            PYBIND11_OVERLOAD(void, Base, save, filename);
        }

        void load(const std::string& filename)
        {
            PYBIND11_OVERLOAD(void, Base, load, filename);
        }
#endif
    };

    template<class Base>
    class PyRandomNumberGeneratorImpl: public PyRandomNumberGenerator<Base>
    {
    public:
        using PyRandomNumberGenerator<Base>::PyRandomNumberGenerator;

        Real random()
        {
            PYBIND11_OVERLOAD(Real, Base, random,);
        }

        Real uniform(Real min, Real max)
        {
            PYBIND11_OVERLOAD(Real, Base, uniform, min, max);
        }

        Integer uniform_int(Integer min, Integer max)
        {
            PYBIND11_OVERLOAD(Integer, Base, uniform_int, min, max);
        }

        Real gaussian(Real sigma, Real mean = 0.0)
        {
            PYBIND11_OVERLOAD(Real, Base, gaussian, sigma, mean);
        }

        Integer binomial(Real p, Integer n)
        {
            PYBIND11_OVERLOAD(Integer, Base, binomial, p, n);
        }

        Real3 direction3d(Real length = 1.0)
        {
            PYBIND11_OVERLOAD(Real3, Base, direction3d, length);
        }

        void seed(Integer val)
        {
            PYBIND11_OVERLOAD(void, Base, seed, val);
        }

        void seed()
        {
            PYBIND11_OVERLOAD(void, Base, seed,);
        }

#ifdef WITH_HDF5
        void save(H5::H5Location* root) const
        {
            PYBIND11_OVERLOAD(void, Base, save, root);
        }

        void load(const H5::H5Location& root)
        {
            PYBIND11_OVERLOAD(void, Base, load, root);
        }

        void save(const std::string& filename) const
        {
            PYBIND11_OVERLOAD(void, Base, save, filename);
        }

        void load(const std::string& filename)
        {
            PYBIND11_OVERLOAD(void, Base, load, filename);
        }
#endif
    };

} // python_api

} // ecell4

#endif /* ECELL4_PYTHON_API_RANDOM_NUMBER_GENERATOR_HPP */
