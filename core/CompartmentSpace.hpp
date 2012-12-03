#ifndef __COMPARTMENT_SPACE_HPP
#define __COMPARTMENT_SPACE_HPP

#include "get_mapper_mf.hpp"
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

    typedef std::vector<Integer> num_molecules_container_type;
    typedef std::vector<Species> species_container_type;

    typedef utils::get_mapper_mf<
        Species, num_molecules_container_type::size_type>::type species_map_type;

    Real volume_;

    num_molecules_container_type num_molecules_;
    species_container_type species_;
    species_map_type index_map_;
};

} // ecell4

#endif /* __COMPARTMENT_SPACE_HPP */
