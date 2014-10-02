#ifndef WORLD_HPP
#define WORLD_HPP


#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/Context.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include "./ParticleTraits.hpp" // This refers ecell4::Particle

#include "ParticleContainerBase.hpp"


#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include "exceptions.hpp"
#include "generator.hpp"
#include "filters.hpp"
//#include "ParticleID.hpp"
#include "SpeciesTypeID.hpp"
#include "SpeciesInfo.hpp"
#include "SerialIDGenerator.hpp"
#include "Transaction.hpp"
#include "Structure.hpp"
#include "Surface.hpp"
#include "Region.hpp"
#include "geometry.hpp"
//#include "GSLRandomNumberGenerator.hpp"
#include "Point.hpp" // XXX: workaround. should be removed later.
#include "utils/pair.hpp"


// For twofold_container
inline
bool is_initialized(std::string const &obj)
{
    return (0 < obj.size());
}

template<typename Tderived_, typename TD_>
struct WorldTraitsBase
{
    // typedef ecell4::Integer size_type;
    typedef std::size_t size_type;
    typedef ecell4::Real length_type;
    typedef ecell4::Real D_type;
    typedef TD_ v_type;
    typedef ecell4::ParticleID particle_id_type;
    typedef ecell4::SerialIDGenerator<particle_id_type> particle_id_generator;
    typedef ecell4::Species::serial_type species_id_type; // std::string
    typedef ecell4::Particle particle_type;
    typedef Sphere particle_shape_type;
    typedef std::string structure_id_type;
    typedef SpeciesInfo<species_id_type, D_type, length_type, structure_id_type> species_type;
    typedef species_type molecule_info;
    typedef ecell4::Position3 point_type;
    typedef typename particle_type::position_type position_type;
    typedef ecell4::GSLRandomNumberGenerator rng_type;
    typedef Structure<Tderived_> structure_type;

    static const Real TOLERANCE = 1e-7;
};

template<typename TD_>
struct WorldTraits: public WorldTraitsBase<WorldTraits<TD_>, TD_>
{
public:
    typedef WorldTraitsBase<WorldTraits<TD_>, TD_> base_type;
    typedef typename base_type::length_type length_type;
    typedef typename base_type::position_type position_type;

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
    static void each_neighbor(Toc_& oc, Tfun_& fun, Tsphere_ const& pos)
    {
        oc.each_neighbor(oc.index(pos), fun);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void each_neighbor(Toc_ const& oc, Tfun_& fun, Tsphere_ const& pos)
    {
        oc.each_neighbor(oc.index(pos), fun);
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

template<typename TD_>
struct CyclicWorldTraits: public WorldTraitsBase<CyclicWorldTraits<TD_>, TD_>
{
public:
    typedef WorldTraitsBase<CyclicWorldTraits<TD_>, TD_> base_type;
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
    static void each_neighbor(Toc_& oc, Tfun_& fun, Tsphere_ const& pos)
    {
        oc.each_neighbor_cyclic(oc.index(pos), fun);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void each_neighbor(Toc_ const& oc, Tfun_& fun, Tsphere_ const& pos)
    {
        oc.each_neighbor_cyclic(oc.index(pos), fun);
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
    typedef ParticleContainer<traits_type> particle_container_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::species_type species_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::particle_id_generator particle_id_generator;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::particle_shape_type particle_shape_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef typename traits_type::rng_type rng_type;
    typedef std::pair<const particle_id_type, particle_type> particle_id_pair;

protected:
    typedef std::map<species_id_type, species_type> species_map;
    typedef std::map<structure_id_type, boost::shared_ptr<structure_type> > structure_map;
    typedef std::set<particle_id_type> particle_id_set;
    typedef std::map<species_id_type, particle_id_set> per_species_particle_id_set;
    typedef select_second<typename species_map::value_type> species_second_selector_type;
    typedef select_second<typename structure_map::value_type> surface_second_selector_type;

public:
    typedef boost::transform_iterator<species_second_selector_type,
            typename species_map::const_iterator> species_iterator;
    typedef boost::transform_iterator<surface_second_selector_type,
            typename structure_map::const_iterator> surface_iterator;
    typedef sized_iterator_range<species_iterator> species_range;
    typedef sized_iterator_range<surface_iterator> structures_range;

protected:

    /** ecell4
     */
    struct MoleculeInfo
    {
        const Real radius;
        const Real D;
        const std::string structure_id;
    };

public:

    World(
        length_type world_size = 1., size_type size = 1)
        : base_type(world_size, size), edge_lengths_(world_size, world_size, world_size)
    {
        rng_ = boost::shared_ptr<rng_type>(new rng_type());
        (*rng_).seed();
    }

    World(
        length_type world_size, size_type size,
        const boost::shared_ptr<rng_type>& rng)
        : base_type(world_size, size), edge_lengths_(world_size, world_size, world_size),
        rng_(rng)
    {
        ;
    }

    virtual particle_id_pair new_particle(species_id_type const& sid,
            position_type const& pos)
    {
        species_type const& species(get_species(sid));
        particle_id_pair retval(pidgen_(),
            particle_type(sid, pos, species.radius(), species.D() ));
        update_particle(retval);
        return retval;
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        typename base_type::particle_matrix_type::iterator i(
                base_type::pmat_.find(pi_pair.first));
        if (i != base_type::pmat_.end())
        {
            if ((*i).second.sid() != pi_pair.second.sid())
            {
                particle_pool_[(*i).second.sid()].erase((*i).first);
                particle_pool_[pi_pair.second.sid()].insert(pi_pair.first);
            }
            base_type::pmat_.update(i, pi_pair);
            return false;
        }
        BOOST_ASSERT(base_type::update_particle(pi_pair));
        particle_pool_[pi_pair.second.sid()].insert(pi_pair.first);
        return true;
    }

    virtual bool remove_particle(particle_id_type const& id)
    {
        bool found(false);
        particle_id_pair pp(this->get_particle(id, found));
        if (!found)
        {
            return false;
        }
        particle_pool_[pp.second.sid()].erase(id);
        base_type::remove_particle(id);
        return true;
    }

    void add_species(species_type const& species)
    {
        species_map_[species.id()] = species;
        particle_pool_[species.id()] = particle_id_set();
    }

    // void add_species(
    //     species_id_type const &sid,
    //     MoleculeInfo const &info,
    //     structure_id_type structure_id = structure_id_type("world"))
    // {
    //     species_type sp(sid, info.D, info.radius, structure_id);
    //     species_map_[sp.id()] = sp;
    //     particle_pool_[sp.id()] = particle_id_set();
    // }

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

    bool add_structure(boost::shared_ptr<structure_type> surface)
    {
        return structure_map_.insert(std::make_pair(surface->id(), surface)).second;
    }

    virtual boost::shared_ptr<structure_type> get_structure(structure_id_type const& id) const
    {
        typename structure_map::const_iterator i(structure_map_.find(id));
        if (structure_map_.end() == i)
        {
            throw not_found(std::string("Unknown surface (id=") + boost::lexical_cast<std::string>(id) + ")");
        }
        return (*i).second;
    }

    structures_range get_structures() const
    {
        return structures_range(
            surface_iterator(structure_map_.begin(), surface_second_selector_type()),
            surface_iterator(structure_map_.end(), surface_second_selector_type()),
            structure_map_.size());
    }

    particle_id_set get_particle_ids(species_id_type const& sid) const
    {
        typename per_species_particle_id_set::const_iterator i(
            particle_pool_.find(sid));
        if (i == particle_pool_.end())
        {
            throw not_found(std::string("Unknown species (id=") + boost::lexical_cast<std::string>(sid) + ")");
        }
        return (*i).second;
    }

    /** ecell4::Space
     */

    inline boost::shared_ptr<rng_type>& rng()
    {
        return rng_;
    }

    // virtual void save(const std::string& filename) const
    // {
    //     throw NotSupported(
    //         "save(const std::string) is not supported by this space class");
    // }

    // virtual void load(const std::string& filename)
    // {
    //     throw NotSupported(
    //         "load(const std::string) is not supported by this space class");
    // }

    virtual const length_type volume() const
    {
        const length_type L(this->world_size());
        return L * L * L;
    }

    virtual const position_type& edge_lengths() const
    {
        return edge_lengths_;
    }

    virtual bool has_particle(const particle_id_type& pid) const
    {
        return base_type::has_particle(pid);
    }

    virtual ecell4::Integer num_particles() const
    {
        return base_type::num_particles();
    }

    virtual ecell4::Integer num_particles_exact(const ecell4::Species& sp) const
    {
        typename per_species_particle_id_set::const_iterator
            i(particle_pool_.find(sp.serial()));
        if (i == particle_pool_.end())
        {
            return 0;
        }
        return (*i).second.size();
    }

    virtual ecell4::Integer num_particles(const ecell4::Species& sp) const
    {
        ecell4::Integer retval(0);
        ecell4::SpeciesExpressionMatcher sexp(sp);
        for (typename per_species_particle_id_set::const_iterator
            i(particle_pool_.begin()); i != particle_pool_.end(); ++i)
        {
            const ecell4::Species tgt((*i).first);
            if (sexp.match(tgt))
            {
                retval += (*i).second.size();
            }
        }
        return retval;
    }

    virtual ecell4::Real get_value(const ecell4::Species& sp) const
    {
        return static_cast<ecell4::Real>(num_molecules(sp));
    }

    virtual ecell4::Real get_value_exact(const ecell4::Species& sp) const
    {
        return static_cast<ecell4::Real>(num_molecules_exact(sp));
    }

    virtual ecell4::Integer num_molecules(const ecell4::Species& sp) const
    {
        ecell4::Integer retval(0);
        ecell4::SpeciesExpressionMatcher sexp(sp);
        for (typename per_species_particle_id_set::const_iterator
            i(particle_pool_.begin()); i != particle_pool_.end(); ++i)
        {
            const ecell4::Species tgt((*i).first);
            retval += sexp.count(tgt) * (*i).second.size();
        }
        return retval;
    }

    virtual ecell4::Integer num_molecules_exact(const ecell4::Species& sp) const
    {
        return num_particles_exact(sp);
    }

    virtual ecell4::Integer num_species() const
    {
        return particle_pool_.size();
    }

    virtual bool has_species(const ecell4::Species& sp) const
    {
        return (particle_pool_.find(sp.serial()) != particle_pool_.end());
    }

    virtual std::vector<std::pair<particle_id_type, particle_type> > list_particles() const
    {
        std::vector<std::pair<particle_id_type, particle_type> > retval;
        retval.reserve(num_particles());
        BOOST_FOREACH(particle_id_pair p, this->get_particles_range())
        {
            retval.push_back(p);
        }
        return retval;
    }

    virtual std::vector<std::pair<particle_id_type, particle_type> >
        list_particles(const ecell4::Species& sp) const
    {
        std::vector<std::pair<particle_id_type, particle_type> > retval;
        ecell4::SpeciesExpressionMatcher sexp(sp);
        for (typename per_species_particle_id_set::const_iterator
            i(particle_pool_.begin()); i != particle_pool_.end(); ++i)
        {
            const ecell4::Species tgt((*i).first);
            if (sexp.match(tgt))
            {
                for (typename particle_id_set::const_iterator j((*i).second.begin());
                    j != (*i).second.end(); ++j)
                {
                    const particle_id_type& pid(*j);
                    retval.push_back(this->get_particle(pid));
                }
            }
        }
        return retval;
    }

    virtual std::vector<std::pair<particle_id_type, particle_type> >
        list_particles_exact(const ecell4::Species& sp) const
    {
        std::vector<std::pair<particle_id_type, particle_type> > retval;
        typename per_species_particle_id_set::const_iterator
            i(particle_pool_.find(sp.serial()));
        if (i == particle_pool_.end())
        {
            return retval;
        }

        for (typename particle_id_set::const_iterator j((*i).second.begin());
            j != (*i).second.end(); ++j)
        {
            const particle_id_type& pid(*j);
            retval.push_back(this->get_particle(pid));
        }
        return retval;
    }

    void bind_to(boost::shared_ptr<ecell4::Model> model)
    {
        if (boost::shared_ptr<ecell4::Model> bound_model = lock_model())
        {
            if (bound_model.get() != model.get())
            {
                std::cerr << "Warning: Model already bound to BDWorld"
                    << std::endl;
            }
        }
        model_ = model;
    }

    boost::shared_ptr<ecell4::Model> lock_model() const
    {
        return model_.lock();
    }

    /**
     * This is a function in the traits of ecell4::ParticleSpace.
     * Be carefull about the difference from
     * "particle_id_pair new_particle(species_id_type const&, position_type const&)".
     */
    std::pair<std::pair<particle_id_type, particle_type>, bool>
    new_particle(const ecell4::Species& sp, const position_type& pos)
    {
        const species_id_type sid(sp.serial());
        species_type const& spinfo(get_species(sid));
        return new_particle(particle_type(sid, pos, spinfo.radius(), spinfo.D()));
    }

    std::pair<std::pair<particle_id_type, particle_type>, bool>
    new_particle(const particle_type& p)
    {
        particle_id_pair retval(pidgen_(), p);
        if (!base_type::check_overlap(
            particle_shape_type(p.position(), p.radius())))
        {
            update_particle(retval);
            return std::make_pair(retval, true);
        }
        else
        {
            return std::make_pair(retval, false);
        }
    }

    /**
     * draw attributes of species and return it as a molecule info.
     * @param sp a species
     * @return info a molecule info
     */
    MoleculeInfo get_molecule_info(const ecell4::Species& sp) const
    {
        Real radius(0.0), D(0.0);
        std::string structure_id("world");

        if (sp.has_attribute("radius") && sp.has_attribute("D"))
        {
            radius = std::atof(sp.get_attribute("radius").c_str());
            D = std::atof(sp.get_attribute("D").c_str());
            if (sp.has_attribute("structure_id"))
            {
                structure_id = sp.get_attribute("structure_id");
            }
        }
        else if (boost::shared_ptr<ecell4::Model> bound_model = lock_model())
        {
            ecell4::Species attributed(bound_model->apply_species_attributes(sp));

            if (attributed.has_attribute("radius")
                && attributed.has_attribute("D"))
            {
                radius = std::atof(
                    attributed.get_attribute("radius").c_str());
                D = std::atof(attributed.get_attribute("D").c_str());
            }

            if (sp.has_attribute("structure_id"))
            {
                structure_id = attributed.get_attribute("structure_id");
            }
        }

        MoleculeInfo info = {radius, D, structure_id};
        return info;
    }

    // species_type get_molecule_info(const ecell4::Species &sp) const
    // {
    //     const ecell4::Species::serial_type sid(sp.serial());
    //     const Real D(std::atof(sp.get_attribute("D").c_str()));
    //     const Real radius(std::atof(sp.get_attribute("radius").c_str()));
    //     const typename species_type::structure_id_type
    //         structure_id(
    //                 sp.has_attribute("structure_id")?
    //                   sp.get_attribute("structure_id"):
    //                   structure_id_type("world"));
    //     return species_type(sid, D, radius, structure_id);
    // }

    // MoleculeInfo get_molecule_info(const ecell4::Species &sp) const
    // {
    //     const Real radius(std::atof(sp.get_attribute("radius").c_str()));
    //     const Real D(std::atof(sp.get_attribute("D").c_str()));
    //     MoleculeInfo info = {radius, D};
    //     return info;
    // }

    /** an adapter function to "void add_species(species_type const&)".
     */
    void add_species(const ecell4::Species& sp)
    {
        MoleculeInfo info(get_molecule_info(sp));
        species_type spinfo(sp.serial(), info.D, info.radius,
            static_cast<typename species_type::structure_id_type>(info.structure_id));
        this->add_species(spinfo);
    }

private:
    particle_id_generator pidgen_;
    species_map species_map_;
    structure_map structure_map_;
    per_species_particle_id_set particle_pool_;

    /** ecell4::Space
     */
    position_type edge_lengths_;
    boost::shared_ptr<rng_type> rng_;
    boost::weak_ptr<ecell4::Model> model_;
};

#endif /* WORLD_HPP */
