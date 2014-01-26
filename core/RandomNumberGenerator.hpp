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

    __attribute__((deprecated)) 
    virtual Real normal(Real mean, Real sigma)  = 0;

    virtual void dir_2d(Real *x, Real *y) = 0;
    virtual void dir_3d(Real *x, Real *y, Real *z) = 0;
    virtual double operator() () = 0;
    virtual void seed(Integer val) = 0;
    virtual void seed() = 0;

    virtual void save(H5::H5File* fout, const std::string& hdf5path) const = 0;
    virtual void load(H5::H5File* fout, const std::string& hdf5path) = 0;
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

    __attribute__((deprecated))
    Real normal(Real loc, Real scale)
    {   // This function is implecated for comatible for epdp::GSLRandomNumberGenerator.
        // This function is the same as uniform().
        std::cout << "hoge" << std::endl;
        return this->uniform(loc, scale);

    }

    Real gaussian(Real mean, Real sigma)
    {
        return gsl_ran_gaussian(rng_.get(), sigma) + mean;
    }

    Integer get_raw()
    {   // epdp
        return gsl_rng_get(rng_.get());
    }

    void dir_2d(Real *x, Real *y)
    {
        gsl_ran_dir_2d(rng_.get(), x, y);
    }

    void dir_3d(Real *x, Real *y, Real *z)
    {
        gsl_ran_dir_3d(rng_.get(), x, y, z);
    }

    Real operator()()
    {
        gsl_rng_uniform(rng_.get());
    }

    void seed(Integer val)
    {   //epdp
        gsl_rng_set(rng_.get(), val);
    }

    void seed()
    {
        gsl_rng_set(rng_.get(), unsigned(std::time(0)));
    }

    void save(H5::H5File* fout, const std::string& hdf5path) const;
    void load(H5::H5File* fout, const std::string& hdf5path);

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
