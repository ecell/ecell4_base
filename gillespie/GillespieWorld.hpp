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
        Real const& volume, boost::shared_ptr<RandomNumberGenerator> rng)
        : cs_(new CompartmentSpaceVectorImpl(volume)), rng_(rng)
    {
        ;
    }

    // SpaceTraits

    Real const& t(void) const;
    void set_t(Real const& t);

    // CompartmentSpaceTraits

    Real const& volume() const
    {
        return cs_->volume();
    }

    Integer num_species(void) const;
    bool has_species(Species const& sp) const;
    Integer num_molecules(Species const& sp) const;

    // CompartmentSpace member functions

    void set_volume(Real const& volume)
    {
        (*cs_).set_volume(volume);
    }

    void add_species(Species const& sp);
    void remove_species(Species const& sp);
    void add_molecules(Species const& sp, Integer const& num);
    void remove_molecules(Species const& sp, Integer const& num);

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
