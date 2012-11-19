#ifndef __RANDOM_NUMBER_GENERATOR_HPP
#define __RANDOM_NUMBER_GENERATOR_HPP

#include <boost/shared_ptr.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "types.hpp"


namespace ecell4
{

class RandomNumberGenerator
{
public:

    virtual Real uniform(Real min, Real max) = 0;
    virtual void seed(Integer val) = 0;
};

class GSLRandomNumberGenerator
    : public RandomNumberGenerator
{
public:

    typedef boost::shared_ptr<gsl_rng> rng_handle;

    Real uniform(Real min, Real max)
    {
        return gsl_rng_uniform(rng_.get()) * (max - min) + min;
    }

    void seed(Integer val)
    {
        gsl_rng_set(rng_.get(), val);
    }

    GSLRandomNumberGenerator(rng_handle hdl)
        : rng_(hdl)
    {
        ;
    }

    GSLRandomNumberGenerator(gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937))
        : rng_(rng, &gsl_rng_free)
    {
        ;
    }

    rng_handle rng_;
};

} // ecell4

#endif /* __RANDOM_NUMBER_GENERATOR_HPP */
