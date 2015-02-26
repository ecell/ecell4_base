#ifndef __ECELL4_RANDOM_NUMBER_GENERATOR_HPP
#define __ECELL4_RANDOM_NUMBER_GENERATOR_HPP

#include <ctime>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef HDF5
#include <hdf5.h>
#include <H5Cpp.h>
#endif

#include "types.hpp"
#include "Real3.hpp"


namespace ecell4
{

class RandomNumberGenerator
{
public:

    virtual ~RandomNumberGenerator()
    {
        ;
    }

    virtual Real random() = 0;
    virtual Real uniform(Real min, Real max) = 0;
    virtual Integer uniform_int(Integer min, Integer max) = 0;
    virtual Real gaussian(Real sigma, Real mean = 0.0) = 0;
    virtual Integer binomial(Real p, Integer n) = 0;
    virtual Real3 direction3d(Real length = 1.0) = 0;

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

    Real random();
    Real uniform(Real min, Real max);
    Integer uniform_int(Integer min, Integer max);
    Real gaussian(Real sigma, Real mean = 0.0);
    Integer binomial(Real p, Integer n);
    Real3 direction3d(Real length);
    void seed(Integer val);
    void seed();

#ifndef HDF5
    void save(H5::CommonFG* root) const;
    void load(const H5::CommonFG& root);
#endif

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

protected:

    rng_handle rng_;
};

} // ecell4

#endif /* __ECELL4_RANDOM_NUMBER_GENERATOR_HPP */
