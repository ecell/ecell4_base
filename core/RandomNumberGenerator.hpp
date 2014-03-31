#ifndef __ECELL4_RANDOM_NUMBER_GENERATOR_HPP
#define __ECELL4_RANDOM_NUMBER_GENERATOR_HPP

#include <ctime>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"


namespace ecell4
{

class RandomNumberGenerator
{
public:

    virtual Real uniform(Real min, Real max) = 0;
    virtual Integer uniform_int(Integer min, Integer max) = 0;
    virtual Real gaussian(Real mean, Real sigma) = 0;

    virtual void seed(Integer val) = 0;
    virtual void seed() = 0;

    virtual void save(H5::CommonFG* root) const = 0;
    virtual void load(const H5::CommonFG& root) = 0;
};

template<typename Telem_>
inline void shuffle(RandomNumberGenerator& rng, std::vector<Telem_>& cont)
{
    typedef std::vector<Telem_> container_type;
    for (typename container_type::size_type i(cont.size()); i > 0;)
    {
        --i;
        typename container_type::size_type const j(rng.uniform_int(0, i));
        std::swap(cont[i], cont[j]);
    }
}

class GSLRandomNumberGenerator
    : public RandomNumberGenerator
{
public:

    typedef boost::shared_ptr<gsl_rng> rng_handle;

public:

    Real uniform(Real min, Real max)
    {
        return gsl_rng_uniform(rng_.get()) * (max - min) + min;
    }

    Integer uniform_int(Integer min, Integer max)
    {
        return gsl_rng_uniform_int(rng_.get(), max - min + 1) + min;
    }

    Real gaussian(Real mean, Real sigma)
    {
        return gsl_ran_gaussian(rng_.get(), sigma) + mean;
    }

    void seed(Integer val)
    {
        gsl_rng_set(rng_.get(), val);
    }

    void seed()
    {
        gsl_rng_set(rng_.get(), unsigned(std::time(0)));
    }

    void save(H5::CommonFG* root) const;
    void load(const H5::CommonFG& root);

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

    inline rng_handle handle()
    {
        return rng_;
    }

protected:

    rng_handle rng_;
};

} // ecell4

#endif /* __ECELL4_RANDOM_NUMBER_GENERATOR_HPP */
