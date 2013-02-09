#ifndef __GILLESPIEWORLD_HPP
#define __GILLESPIEWORLD_HPP

#include <stdexcept>
#include <map>
#include <boost/scoped_ptr.hpp>
#include <string>

#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/CompartmentSpace.hpp>
#include <ecell4/core/Species.hpp>


namespace ecell4
{

namespace gillespie
{

class GillespieWorld
{
public:

    GillespieWorld(
        const Real& volume, boost::shared_ptr<RandomNumberGenerator> rng)
        : cs_(new CompartmentSpaceVectorImpl(volume)), rng_(rng)
    {
        ;
    }

    // SpaceTraits

    const Real& t(void) const;
    void set_t(const Real& t);

    // CompartmentSpaceTraits

    const Real& volume() const
    {
        return cs_->volume();
    }

    Integer num_species(void) const;
    bool has_species(const Species& sp) const;
    Integer num_molecules(const Species& sp) const;

    // CompartmentSpace member functions

    void set_volume(const Real& volume)
    {
        (*cs_).set_volume(volume);
    }

    void add_species(const Species& sp);
    void remove_species(const Species& sp);
    void add_molecules(const Species& sp, const Integer& num);
    void remove_molecules(const Species& sp, const Integer& num);

    // Optional members

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return rng_;
    }

private:

    boost::scoped_ptr<CompartmentSpace> cs_;
    boost::shared_ptr<RandomNumberGenerator> rng_;
};

} // gillespie

} // ecell4


#endif // __GILLESPIEWORLD_HPP
