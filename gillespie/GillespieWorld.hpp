#ifndef __GILLESPIEWORLD_HPP
#define __GILLESPIEWORLD_HPP

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
        Real const &volume, boost::shared_ptr<RandomNumberGenerator> rng)
        : cs_(new CompartmentSpaceVectorImpl(volume)), rng_(rng)
    {
        ;
    }
    // about time
    void set_t(Real const &t);
    Real t(void);

    Real const& volume() const
    {
        return cs_->volume();
    }

    // about molecules states
    // immutable functions.
    Integer num_species(void);
    bool has_species(Species const &sp);
    Integer num_molecules(Species const& sp);

    // mutable functions.
    void add_species(Species const &sp);
    void remove_species(Species const &sp);
    void add_molecules(Species const &sp, Integer const &num);
    // I think it is better that the name of this function is 'decrease_molecules()'.
    void remove_molecules(Species const &sp, Integer const &num);

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
