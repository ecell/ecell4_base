#ifndef GSLRANDOMNUMBERGENERATOR_HPP
#define GSLRANDOMNUMBERGENERATOR_HPP

#include <boost/shared_ptr.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class GSLRandomNumberGenerator
{
public:
    typedef boost::shared_ptr<gsl_rng> rng_handle;

    double normal(double loc, double scale)
    {
        return gsl_ran_gaussian(rng_.get(), scale) + loc;
    }

    unsigned long int get_raw()
    {
        return gsl_rng_get(rng_.get());
    }

    double uniform(double min, double max)
    {
        return gsl_rng_uniform(rng_.get()) * (max - min) + min;
    }

    int uniform_int(int min, int max)
    {
        return gsl_rng_uniform_int(rng_.get(), max - min + 1) + min;
    }

    void dir_2d(double *x, double *y)
    {
        gsl_ran_dir_2d(rng_.get(), x, y);
    }

    void dir_3d(double *x, double *y, double *z)
    {
        gsl_ran_dir_3d(rng_.get(), x, y, z);
    }

    double operator()()
    {
        return gsl_rng_uniform(rng_.get());
    }

    void seed(unsigned long int val)
    {
        gsl_rng_set(rng_.get(), val);
    }

    GSLRandomNumberGenerator(rng_handle hdl): rng_(hdl)
    {
    }

    GSLRandomNumberGenerator(gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937)): rng_(rng, &gsl_rng_free) {}

    rng_handle rng_;
};

#endif /* GSLRANDOMNUMBERGENERATOR_HPP */
