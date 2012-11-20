#ifndef __RANDOM_NUMBER_GENERATOR_HPP
#define __RANDOM_NUMBER_GENERATOR_HPP

#include <vector>
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
    virtual Integer uniform_int(Integer min, Integer max) = 0;
    virtual void seed(Integer val) = 0;
};

template <typename Telem_>
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

    Real uniform(Real min, Real max)
    {
        return gsl_rng_uniform(rng_.get()) * (max - min) + min;
    }

    Integer uniform_int(Integer min, Integer max)
    {
        return gsl_rng_uniform_int(rng_.get(), max - min + 1) + min;
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
