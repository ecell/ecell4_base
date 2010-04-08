#ifndef WORLD_HPP
#define WORLD_HPP

#include "ParticleContainerBase.hpp"

#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include "exceptions.hpp"
#include "generator.hpp"
#include "filters.hpp"
#include "Particle.hpp"
#include "ParticleID.hpp"
#include "SpeciesTypeID.hpp"
#include "SpeciesInfo.hpp"
#include "SerialIDGenerator.hpp"
#include "Transaction.hpp"
#include "Structure.hpp"
#include "geometry.hpp"

template<typename Tlen_, typename TD_>
struct WorldTraitsBase
{
    typedef std::size_t size_type;
    typedef Tlen_ length_type;
    typedef TD_ D_type;
    typedef ParticleID particle_id_type;
    typedef SerialIDGenerator<particle_id_type> particle_id_generator;
    typedef SpeciesTypeID species_id_type;
    typedef Particle<length_type, D_type, species_id_type> particle_type;
    typedef std::string surface_id_type;
    typedef SpeciesInfo<species_id_type, D_type, length_type, surface_id_type> species_type;
    typedef Vector3<length_type> point_type;
    typedef typename particle_type::shape_type::position_type position_type;
    typedef Structure<surface_id_type> surface_type;

    template<typename Tval_>
    static Tval_ apply_boundary(Tval_ const& v, length_type const& world_size)
    {
        return v;
    }

    template<typename Tval_>
    static Tval_ cyclic_transpose(Tval_ const& p0, Tval_ const& p1, length_type const& world_size)
    {
        return p0;
    }

    template<typename T1_, typename T2_>
    static length_type distance(T1_ const& p0, T2_ const& p1, length_type const& world_size)
    {
        return ::distance(p0, p1);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void take_neighbor(Toc_& oc, Tfun_& fun, const Tsphere_& cmp)
    {
        take_neighbor(oc, fun, cmp);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void take_neighbor(Toc_ const& oc, Tfun_& fun, const Tsphere_& cmp)
    {
        take_neighbor(oc, fun, cmp);
    }
};

template<typename Tlen_, typename TD_>
struct CyclicWorldTraits: public WorldTraitsBase<Tlen_, TD_>
{
protected:
    typedef WorldTraitsBase<Tlen_, TD_> base_type;
public:
    typedef typename base_type::length_type length_type;
    typedef typename base_type::position_type position_type;

    template<typename Tval_>
    static Tval_ apply_boundary(Tval_ const& v, length_type const& world_size)
    {
        return ::apply_boundary(v, world_size);
    }

    static length_type cyclic_transpose(length_type const& p0, length_type const& p1, length_type const& world_size)
    {
        return ::cyclic_transpose(p0, p1, world_size);
    }

    static position_type cyclic_transpose(position_type const& p0, position_type const& p1, length_type const& world_size)
    {
        return ::cyclic_transpose(p0, p1, world_size);
    }

    template<typename T1_, typename T2_>
    static length_type distance(T1_ const& p0, T2_ const& p1, length_type const& world_size)
    {
        return distance_cyclic(p0, p1, world_size);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void take_neighbor(Toc_& oc, Tfun_& fun, const Tsphere_& cmp)
    {
        take_neighbor_cyclic(oc, fun, cmp);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void take_neighbor(Toc_ const& oc, Tfun_& fun, const Tsphere_& cmp)
    {
        take_neighbor_cyclic(oc, fun, cmp);
    }
};

template<typename Ttraits_>
class World: public ParticleContainerBase<World<Ttraits_>, Ttraits_>
{
public:
    typedef Ttraits_ traits_type;
    typedef ParticleContainerBase<World> base_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::species_type species_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::particle_id_generator particle_id_generator;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::particle_type::shape_type particle_shape_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::surface_id_type surface_id_type;
    typedef typename traits_type::surface_type surface_type;
    typedef std::pair<const particle_id_type, particle_type> particle_id_pair;

protected:
    typedef std::map<species_id_type, species_type> species_map;
    typedef std::map<surface_id_type, surface_type*> surface_map;
    typedef select_second<typename species_map::value_type> species_second_selector_type;
    typedef select_second<typename surface_map::value_type> surface_second_selector_type;

public:
    typedef boost::transform_iterator<species_second_selector_type,
            typename species_map::const_iterator> species_iterator;
    typedef boost::transform_iterator<surface_second_selector_type,
            typename surface_map::const_iterator> surface_iterator;
    typedef sized_iterator_range<species_iterator> species_range;
    typedef sized_iterator_range<surface_iterator> surfaces_range;

public:
    World(length_type world_size = 1., size_type size = 1)
        : base_type(world_size, size) {}

    particle_id_pair new_particle(species_id_type const& sid,
            position_type const& pos)
    {
        species_type const& species(get_species(sid));
        particle_id_pair retval(pidgen_(),
            particle_type(sid, particle_shape_type(pos, species.radius()),
                          species.D()));
        base_type::pmat_.update(retval);
        return retval;
    }

    void add_species(species_type const& species)
    {
        species_map_[species.id()] = species;
    }

    virtual species_type const& get_species(species_id_type const& id) const
    {
        typename species_map::const_iterator i(species_map_.find(id));
        if (species_map_.end() == i)
        {
            throw not_found(std::string("Unknown species (id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return (*i).second;
    }

    species_range get_species() const
    {
        return species_range(
            species_iterator(species_map_.begin(), species_second_selector_type()),
            species_iterator(species_map_.end(), species_second_selector_type()),
            species_map_.size());
    }

    bool add_surface(surface_type* surface)
    {
        return surface_map_.insert(std::make_pair(surface->id(), surface)).second;
    }

    virtual surface_type const& get_surface(surface_id_type const& id) const
    {
        typename surface_map::const_iterator i(surface_map_.find(id));
        if (surface_map_.end() == i)
        {
            throw not_found(std::string("Unknown surface (id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return *(*i).second;
    }

    surfaces_range get_surfaces() const
    {
        return surfaces_range(
            surface_iterator(surface_map_.begin(), surface_second_selector_type()),
            surface_iterator(surface_map_.end(), surface_second_selector_type()),
            surface_map_.size());
    }

private:
    particle_id_generator pidgen_;
    species_map species_map_;
    surface_map surface_map_;
};

#endif /* WORLD_HPP */
