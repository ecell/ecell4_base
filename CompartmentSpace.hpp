#ifndef __COMPARTMENT_SPACE_HPP
#define __COMPARTMENT_SPACE_HPP

#include <map>

#include "types.hpp"
#include "Species.hpp"
#include "Space.hpp"


namespace ecell4
{

class CompartmentSpace
    : public Space
{
public:

    virtual Real const& volume() const = 0;
    virtual void set_volume(Real volume) = 0;

    virtual void add_species(Species const& sp) = 0;
    virtual void remove_species(Species const& sp) = 0;
    virtual Integer num_species() const = 0;

    virtual Integer num_molecules(Species const& sp) const = 0;
    virtual void add_molecules(Species const& sp, Integer const& num) = 0;
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
