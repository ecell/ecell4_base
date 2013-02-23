#ifndef __ECELL4_EGFRD_EGFRD_WORLD_HPP
#define __ECELL4_EGFRD_EGFRD_WORLD_HPP

#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>

// epdp
#include "World.hpp"
#include "ParticleModel.hpp"
#include "SpeciesType.hpp"
#include "SpeciesTypeID.hpp"
#include "CuboidalRegion.hpp"
#include "NetworkRules.hpp"
#include "ReactionRule.hpp"
#include "EGFRDSimulator.hpp"
// epdp

#include <ecell4/core/types.hpp>
#include <ecell4/core/get_mapper_mf.hpp>
#include <ecell4/core/extras.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/Particle.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>


namespace ecell4
{

namespace egfrd
{

struct ParticleInfo
{
    const Real radius;
    const Real D;
};

class EGFRDWorld
{
public:

    typedef ParticleInfo particle_info_type;

    typedef ::World< ::CyclicWorldTraits<Real, Real> > world_type;
    typedef ::ParticleModel particle_model_type;
    typedef EGFRDSimulator< ::EGFRDSimulatorTraitsBase<EGFRDWorld::world_type> >
    simulator_type;

protected:

    typedef std::vector<std::pair<Species::serial_type, ::SpeciesTypeID> >
    species_type_id_container_type;
    typedef ::CuboidalRegion<simulator_type::traits_type> cuboidal_region_type;
    typedef typename world_type::traits_type::structure_id_type
    structure_id_type;
    typedef simulator_type::traits_type::network_rules_type network_rules_type;

public:

    EGFRDWorld(
        const Real& world_size, const Integer& matrix_size,
        boost::shared_ptr<GSLRandomNumberGenerator> rng)
        : world_(new world_type(world_size, matrix_size)),
          rng_(rng), internal_rng_(rng->handle()), t_(0.0)
    {
        const world_type::position_type x(
            translate(divide(edge_lengths(), 2)));
        (*world_).add_structure(
            boost::shared_ptr<cuboidal_region_type>(
                new cuboidal_region_type(
                    "world", typename cuboidal_region_type::shape_type(x, x))));

        ;
    }

    /**
     * create and add a new particle
     * @param p a particle
     * @return pid a particle id
     */
    ParticleID new_particle(const Particle& p)
    {
        // add_species(p.species());

        world_type::particle_id_pair retval(
            (*world_).new_particle(
                find(p.species()), translate(p.position())));
        return translate(retval.first);
    }

    /**
     * draw attributes of species and return it as a particle info.
     * @param sp a species
     * @return info a particle info
     */
    ParticleInfo get_particle_info(const Species& sp) const
    {
        const Real radius(std::atof(sp.get_attribute("radius").c_str()));
        const Real D(std::atof(sp.get_attribute("D").c_str()));
        ParticleInfo info = {radius, D};
        return info;
    }

    // SpaceTraits

    const Real& t() const
    {
        return t_;
    }

    void set_t(const Real& t)
    {
        t_ = t;
    }

    // ParticleSpaceTraits

    Position3 edge_lengths() const
    {
        const double& world_size((*world_).world_size());
        return Position3(world_size, world_size, world_size);
    }

    Integer num_particles() const
    {
        return static_cast<Integer>((*world_).num_particles());
    }

    Integer num_particles(const Species& sp) const
    {
        return static_cast<Integer>(
            (*world_).get_particle_ids(find(sp)).size());
    }

    std::vector<std::pair<ParticleID, Particle> > list_particles() const
    {
        std::vector<std::pair<ParticleID, Particle> > particles(num_particles());
        BOOST_FOREACH(const world_type::particle_id_pair& pid_pair,
                      (*world_).get_particles_range())
        {
            particles.push_back(
                std::make_pair(translate(pid_pair.first),
                               translate(pid_pair.second)));
        }
        return particles;
    }

    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& species) const
    {
        const ::SpeciesTypeID target(find(species));
        // std::vector<std::pair<ParticleID, Particle> > particles(num_particles());
        std::vector<std::pair<ParticleID, Particle> > particles;
        BOOST_FOREACH(const world_type::particle_id_pair& pid_pair,
                      (*world_).get_particles_range())
        {
            if (pid_pair.second.sid() == target)
            {
                particles.push_back(
                    std::make_pair(translate(pid_pair.first),
                                   translate(pid_pair.second)));
            }
        }
        return particles;
    }

    bool has_particle(const ParticleID& pid) const
    {
        return (*world_).has_particle(translate(pid));
    }

    // ParticleSpace member functions

    bool update_particle(const ParticleID& pid, const Particle& p)
    {
        return (*world_).update_particle(
            std::make_pair(translate(pid), translate(p)));
    }

    std::pair<ParticleID, Particle>
    get_particle(const ParticleID& pid) const
    {
        world_type::particle_id_pair
            pid_pair((*world_).get_particle(translate(pid)));
        return std::make_pair(pid, translate(pid_pair.second));
    }

    void remove_particle(const ParticleID& pid)
    {
        (*world_).remove_particle(translate(pid));
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius) const
    {
        boost::scoped_ptr<world_type::particle_id_pair_and_distance_list>
            overlapped(
                (*world_).check_overlap(
                    world_type::particle_shape_type(translate(pos), radius)));
        if (!overlapped || ::size(*overlapped) == 0)
        {
            return std::vector<
                std::pair<std::pair<ParticleID, Particle>, Real> >();
        }

        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
            retval(::size(*overlapped));
        for (world_type::particle_id_pair_and_distance_list::const_iterator
                 i(overlapped->begin()); i != overlapped->end(); ++i)
        {
            retval.push_back(
                std::make_pair(
                    std::make_pair(
                        translate((*i).first.first), translate((*i).first.second)),
                    (*i).second));
        }
        return retval;
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius, const ParticleID& ignore) const
    {
        boost::scoped_ptr<world_type::particle_id_pair_and_distance_list>
            overlapped(
                (*world_).check_overlap(
                    world_type::particle_shape_type(translate(pos), radius),
                    translate(ignore)));
        if (!overlapped || ::size(*overlapped) == 0)
        {
            return std::vector<
                std::pair<std::pair<ParticleID, Particle>, Real> >();
        }

        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
            retval(::size(*overlapped));
        for (world_type::particle_id_pair_and_distance_list::const_iterator
                 i(overlapped->begin()); i != overlapped->end(); ++i)
        {
            retval.push_back(
                std::make_pair(
                    std::make_pair(
                        translate((*i).first.first), translate((*i).first.second)),
                    (*i).second));
        }
        return retval;
    }

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Position3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
    {
        boost::scoped_ptr<world_type::particle_id_pair_and_distance_list>
            overlapped(
                (*world_).check_overlap(
                    world_type::particle_shape_type(translate(pos), radius),
                    translate(ignore1), translate(ignore2)));
        if (!overlapped || ::size(*overlapped) == 0)
        {
            return std::vector<
                std::pair<std::pair<ParticleID, Particle>, Real> >();
        }

        std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
            retval(::size(*overlapped));
        for (world_type::particle_id_pair_and_distance_list::const_iterator
                 i(overlapped->begin()); i != overlapped->end(); ++i)
        {
            retval.push_back(
                std::make_pair(
                    std::make_pair(
                        translate((*i).first.first), translate((*i).first.second)),
                    (*i).second));
        }
        return retval;
    }

    inline Position3 periodic_transpose(
        const Position3& pos1, const Position3& pos2) const
    {
        return translate(
            (*world_).cyclic_transpose(translate(pos1), translate(pos2)));
    }

    inline Position3 apply_boundary(const Position3& pos) const
    {
        return translate((*world_).apply_boundary(translate(pos)));
    }

    inline Real distance_sq(const Position3& pos1, const Position3& pos2) const
    {
        const Real L(distance(pos1, pos2));
        return L * L;
    }

    inline Real distance(const Position3& pos1, const Position3& pos2) const
    {
        return static_cast<Real>(
            (*world_).distance(translate(pos1), translate(pos2)));
    }

    // CompartmentSpaceTraits

    Real volume() const
    {
        const Position3 lengths(edge_lengths());
        return lengths[0] * lengths[1] * lengths[2];
    }

    Integer num_species() const
    {
        return sid_container_.size();
    }

    bool has_species(const Species& sp) const
    {
        return (find_species_type_id(sp) != sid_container_.end());
    }

    Integer num_molecules(const Species& sp) const
    {
        return num_particles(sp);
    }

    // CompartmentSpace member functions

    void add_species(const Species& sp)
    {
        if (!has_species(sp))
        {
            // add ::SpeciesType to ::ParticleModel
            boost::shared_ptr< ::SpeciesType> st(new ::SpeciesType());
            (*st)["name"] = boost::lexical_cast<std::string>(sp.name());
            (*st)["D"] = boost::lexical_cast<std::string>(sp.get_attribute("D"));
            (*st)["radius"] = boost::lexical_cast<std::string>(
                sp.get_attribute("radius"));
            model_.add_species_type(st);

            // create a map between Species and ::SpeciesType
            sid_container_.push_back(std::make_pair(sp.serial(), st->id()));

            // add ::SpeciesInfo to ::World
            const std::string& structure_id((*st)["structure"]);
            (*world_).add_species(
                typename world_type::traits_type::species_type(
                    st->id(),
                    boost::lexical_cast<typename world_type::traits_type::D_type>(
                        (*st)["D"]),
                    boost::lexical_cast<typename world_type::length_type>(
                        (*st)["radius"]),
                    boost::lexical_cast<structure_id_type>(
                        structure_id.empty() ? "world": structure_id)));
        }
    }

    void remove_species(const Species& sp)
    {
        throw NotImplemented("Not implemented yet.");
    }

    void add_molecules(const Species& sp, const Integer& num)
    {
        extras::throw_in_particles(*this, sp, num, *rng());
    }

    // void remove_molecules(const Species& sp, const Integer& num)
    // {
    //     throw NotSupported("Not supported. Use remove_particle instead.");
    // }

    // Optional members

    void add_reaction_rule(const ReactionRule& rr)
    {
        model_.network_rules().add_reaction_rule(translate(rr));
    }

    // const particle_container_type& particles() const
    // {
    //     throw NotSupported("Not supported. Use list_particles() instead.");
    // }

    inline boost::shared_ptr<GSLRandomNumberGenerator> rng()
    {
        return rng_;
    }

    simulator_type* create_simulator(const Integer& dissociation_retry_moves)
    {
        return new simulator_type(
            world_, boost::shared_ptr<network_rules_type>(
                new network_rules_type(model_.network_rules())),
            internal_rng_, dissociation_retry_moves);
    }

protected:

    template<typename Tfirst_, typename Tsecond_>
    struct pair_first_element_unary_predicator
    {
        typedef std::pair<Tfirst_, Tsecond_> element_type;

        pair_first_element_unary_predicator(const Tfirst_& target)
            : target_(target)
        {
            ; // do nothing
        }

        bool operator()(const element_type& v)
        {
            return v.first == target_;
        }

    protected:

        Tfirst_ target_;
    };

    template<typename Tfirst_, typename Tsecond_>
    struct pair_second_element_unary_predicator
    {
        typedef std::pair<Tfirst_, Tsecond_> element_type;

        pair_second_element_unary_predicator(const Tsecond_& target)
            : target_(target)
        {
            ; // do nothing
        }

        bool operator()(const element_type& v)
        {
            return v.second == target_;
        }

    protected:

        Tsecond_ target_;
    };

    species_type_id_container_type::const_iterator
    find_species_type_id(const Species& sp) const
    {
        pair_first_element_unary_predicator<
            Species::serial_type, ::SpeciesTypeID> predicator(sp.serial());
        species_type_id_container_type::const_iterator
            i(std::find_if(
                  sid_container_.begin(), sid_container_.end(), predicator));
        return i;
    }

    inline ::SpeciesTypeID find(const Species& sp) const
    {
        species_type_id_container_type::const_iterator
            i(find_species_type_id(sp));
        if (i == sid_container_.end())
        {
            throw NotFound("Species not found");
        }
        return (*i).second;
    }

    inline Species find(const ::SpeciesTypeID& sid) const
    {
        pair_second_element_unary_predicator<
            Species::serial_type, ::SpeciesTypeID> predicator(sid);
        species_type_id_container_type::const_iterator
            i(std::find_if(
                  sid_container_.begin(), sid_container_.end(), predicator));
        if (i == sid_container_.end())
        {
            throw NotFound("Species not found");
        }
        return (*i).first;
    }

    inline ParticleID translate(const world_type::particle_id_type& pid) const
    {
        return ParticleID(ParticleID::value_type(pid.lot(), pid.serial()));
    }

    inline world_type::particle_id_type translate(const ParticleID& pid) const
    {
        return world_type::particle_id_type(
            world_type::particle_id_type::value_type(pid.lot(), pid.serial()));
    }

    inline world_type::position_type translate(const Position3& pos) const
    {
        return world_type::position_type(pos[0], pos[1], pos[2]);
    }

    inline Position3 translate(const world_type::position_type& pos) const
    {
        return Position3(pos[0], pos[1], pos[2]);
    }

    inline world_type::particle_type translate(const Particle& p) const
    {
        const Real radius(p.radius()), D(p.D());
        world_type::particle_type retval(
            find(p.species()), world_type::particle_shape_type(
                translate(p.position()), radius), D);
        return retval;
    }

    inline Particle translate(const world_type::particle_type& p) const
    {
        const Real radius(p.radius()), D(p.D());
        Particle retval(
            find(p.sid()), translate(p.position()), radius, D);
        return retval;
    }

    inline ReactionRule translate(const ::ReactionRule& rr) const
    {
        throw NotImplemented("Not implemented yet.");
    }

    inline ::ReactionRule translate(const ReactionRule& rr) const
    {
        ::ReactionRule retval;
        std::vector< ::SpeciesTypeID> products;
        for (ReactionRule::product_container_type::const_iterator
                 j(rr.products().begin()); j != rr.products().end(); ++j)
        {
            products.push_back(find(*j));
        }

        ReactionRule::reactant_container_type::const_iterator
            r(rr.reactants().begin());
        switch (rr.reactants().size())
        {
        case 1:
            {
                const ::SpeciesTypeID sid1(find(*r));
                retval = ::new_reaction_rule(sid1, products, rr.k());
            }
            break;
        case 2:
            {
                const ::SpeciesTypeID sid1(find(*r));
                ++r;
                const ::SpeciesTypeID sid2(find(*r));
                retval = ::new_reaction_rule(sid1, sid2, products, rr.k());
            }
            break;
        default:
            throw NotSupported("the number of reactants must be 1 or 2.");
        }

        return retval;
    }

protected:

    boost::shared_ptr<world_type> world_;
    boost::shared_ptr<GSLRandomNumberGenerator> rng_;
    world_type::traits_type::rng_type internal_rng_;
    Real t_;

    particle_model_type model_;
    species_type_id_container_type sid_container_;
};

} // egfrd

} // ecell4

#endif // __ECELL4_EGFRD_EGFRD_WORLD_HPP
