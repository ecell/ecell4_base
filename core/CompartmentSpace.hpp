#ifndef __COMPARTMENT_SPACE_HPP
#define __COMPARTMENT_SPACE_HPP

#include <map>

#include "types.hpp"
#include "exceptions.hpp"
#include "Species.hpp"
#include "Space.hpp"


namespace ecell4
{

class CompartmentSpace
    : public Space
{
public:

    virtual Real const& volume() const
    {
        throw NotImplemented("volume() not implemented");
    }

    virtual Integer num_species() const
    {
        throw NotImplemented("num_species() not implemented");
    }

    virtual bool has_species(Species const& sp) const
    {
        throw NotImplemented("has_species() not implemented");
    }

    virtual Integer num_molecules(Species const& sp) const
    {
        throw NotImplemented("num_molecules() not implemented");
    }

    /**
     * set volume.
     * this function is a member of CompartmentSpace.
     * @param volume a nonzero positive Real value
     */
    virtual void set_volume(Real volume) = 0;

    /**
     * add a species.
     * this function is a member of CompartmentSpace.
     * @param sp a species
     */
    virtual void add_species(Species const& sp) = 0;

    /**
     * remove a species.
     * this function is a member of CompartmentSpace.
     * @param sp a species
     */
    virtual void remove_species(Species const& sp) = 0;

    /**
     * increase the number of molecules.
     * this function is a member of CompartmentSpace.
     * @param sp a species
     * @param num a number of molecules
     */
    virtual void add_molecules(Species const& sp, Integer const& num) = 0;

    /**
     * decrease the number of molecules.
     * this function is a member of CompartmentSpace.
     * @param sp a species
     * @param num a number of molecules
     */
    virtual void remove_molecules(Species const& sp, Integer const& num) = 0;
};

class CompartmentSpaceVectorImpl
    : public CompartmentSpace
{
public:

    typedef std::vector<Integer>::size_type index_type;
    typedef std::map<Species, index_type> index_map_type;

    CompartmentSpaceVectorImpl(Real const& volume)
        : volume_(1)
    {
        set_volume(volume);
    }

    Real const& volume() const;
    void set_volume(Real volume);

    void add_species(Species const& sp);
    void remove_species(Species const& sp);
    Integer num_species() const;
    bool has_species(Species const& sp) const;

    Integer num_molecules(Species const& sp) const;
    void add_molecules(Species const& sp, Integer const& num);
    void remove_molecules(Species const& sp, Integer const& num);

protected:

    Real volume_;

    std::vector<Integer> num_molecules_;
    std::vector<Species> species_;
    index_map_type index_map_;
};

} // ecell4

#endif /* __COMPARTMENT_SPACE_HPP */
