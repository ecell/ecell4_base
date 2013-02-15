#ifndef __SPATIOCYTE_WORLD_HPP
#define __SPATIOCYTE_WORLD_HPP

#include <sstream>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#ifndef __ECELL3_EXCEPTIONS_HPP
#define __ECELL3_EXCEPTIONS_HPP
#undef __EXCEPTIONS_HPP
#include <libecs/Exceptions.hpp>
#endif

#include <libecs/libecs.hpp>
#include <libecs/Entity.hpp>
#include <libecs/Variable.hpp>
#include <libecs/Process.hpp>

#ifndef __ECELL3_MODEL_HPP
#define __ECELL3_MODEL_HPP
#undef __MODEL_HPP
#include <libecs/Model.hpp>
#endif

#include <SpatiocyteCommon.hpp>

#include <ecell4/core/extras.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/ParticleSpace.hpp>


namespace ecell4
{

namespace spatiocyte
{

struct ParticleInfo
{
    const Real radius;
    const Real D;
};

class SpatiocyteWorld
{
public:

    typedef ParticleInfo particle_info_type;

protected:

    typedef std::vector<std::pair<Species::serial_type, libecs::String> >
    species_container_type;

public:

    SpatiocyteWorld(const Position3& edge_lengths, const Real& voxel_radius)
        : edge_lengths_(edge_lengths), voxel_radius_(voxel_radius), t_(0.0),
          num_reactions_(0)
    {
        libecs::initialize();
        model_ = new libecs::Model(*libecs::createDefaultModuleMaker());
        setup_model();
    }

    ~SpatiocyteWorld()
    {
        delete model_;
        libecs::finalize();
    }

    // /**
    //  * create and add a new particle
    //  * @param p a particle
    //  * @return pid a particle id
    //  */
    // ParticleID new_particle(const Particle& p)
    // {
    //     ParticleID pid(pidgen_());
    //     (*ps_).update_particle(pid, p);
    //     return pid;
    // }

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
        if (t < 0.0)
        {
            throw std::invalid_argument("the time must be positive.");
        }
        t_ = t;
    }

    // ParticleSpaceTraits

    const Position3& edge_lengths() const
    {
        return edge_lengths_;
    }

    // Integer num_particles() const
    // {
    //     return (*ps_).num_particles();
    // }

    Integer num_particles(const Species& sp) const
    {
        species_container_type::const_iterator i(find_variable_fullid(sp));
        if (i == species_.end())
        {
            return 0;
        }

        return static_cast<Integer>(get_variable((*i).second)->getValue());
    }

    // bool has_particle(const ParticleID& pid) const
    // {
    //     return (*ps_).has_particle(pid);
    // }

    // std::vector<std::pair<ParticleID, Particle> > list_particles() const
    // {
    //     return (*ps_).list_particles();
    // }

    // std::vector<std::pair<ParticleID, Particle> >
    // list_particles(const Species& species) const
    // {
    //     return (*ps_).list_particles(species);
    // }

    // // ParticleSpace member functions

    // bool update_particle(const ParticleID& pid, const Particle& p)
    // {
    //     return (*ps_).update_particle(pid, p);
    // }

    // std::pair<ParticleID, Particle>
    // get_particle(const ParticleID& pid) const
    // {
    //     return (*ps_).get_particle(pid);
    // }

    // void remove_particle(const ParticleID& pid)
    // {
    //     (*ps_).remove_particle(pid);
    // }

    // std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    // list_particles_within_radius(
    //     const Position3& pos, const Real& radius) const
    // {
    //     return (*ps_).list_particles_within_radius(pos, radius);
    // }

    // std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    // list_particles_within_radius(
    //     const Position3& pos, const Real& radius, const ParticleID& ignore) const
    // {
    //     return (*ps_).list_particles_within_radius(pos, radius, ignore);
    // }

    // std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    // list_particles_within_radius(
    //     const Position3& pos, const Real& radius,
    //     const ParticleID& ignore1, const ParticleID& ignore2) const
    // {
    //     return (*ps_).list_particles_within_radius(pos, radius, ignore1, ignore2);
    // }

    // inline Position3 periodic_transpose(
    //     const Position3& pos1, const Position3& pos2) const
    // {
    //     return (*ps_).periodic_transpose(pos1, pos2);
    // }

    // inline Position3 apply_boundary(const Position3& pos) const
    // {
    //     return (*ps_).apply_boundary(pos);
    // }

    // inline Real distance_sq(const Position3& pos1, const Position3& pos2) const
    // {
    //     return (*ps_).distance_sq(pos1, pos2);
    // }

    // inline Real distance(const Position3& pos1, const Position3& pos2) const
    // {
    //     return (*ps_).distance(pos1, pos2);
    // }

    // CompartmentSpaceTraits

    Real volume() const
    {
        const Position3 lengths(edge_lengths());
        return lengths[0] * lengths[1] * lengths[2];
    }

    Integer num_molecules(const Species& sp) const
    {
        return num_particles(sp);
    }

    bool has_species(const Species& sp) const
    {
        return (find_variable_fullid(sp) != species_.end());
    }

    void add_molecules(const Species& sp, const Integer& num)
    {
        species_container_type::const_iterator i(find_variable_fullid(sp));
        if (i == species_.end())
        {
            add_species(sp);
            i = find_variable_fullid(sp);
        }

        libecs::Variable* variable_ptr(get_variable((*i).second));
        variable_ptr->setValue(variable_ptr->getValue() + num);
    }

    // CompartmentSpace member functions

    void add_species(const Species& sp)
    {
        if (!has_species(sp))
        {
            libecs::String fullid_str("Variable:/:" + sp.name());
            species_.push_back(std::make_pair(sp.serial(), fullid_str));

            create_variable(fullid_str, 0);
            const ParticleInfo info(get_particle_info(sp));
            create_diffusion_process(fullid_str, info.D);

            get_molecule_populate_process()->registerVariableReference(
                "_", fullid_str, 0);
        }
    }

    // Optional members

    void step()
    {
        (*model_).step();
    }

    Real dt() const
    {
        return static_cast<Real>((*model_).getStepper("SS")->getStepInterval());
    }

    void add_reaction_rule(const ReactionRule& rr)
    {
        switch (rr.reactants().size())
        {
        case 1:
            {
                std::stringstream ss;
                ss << "Process:/:rr" << ++num_reactions_;
                libecs::Process* const process_ptr(
                    dynamic_cast<libecs::Process*>(
                        (*model_).createEntity(
                            "SpatiocyteNextReactionProcess",
                            libecs::FullID(ss.str()))));
                process_ptr->loadProperty("k", libecs::Polymorph(rr.k()));

                ReactionRule::reactant_container_type::const_iterator
                    r(rr.reactants().begin());
                process_ptr->registerVariableReference(
                    "_", get_variable_fullid(*r), -1);
                for (ReactionRule::reactant_container_type::const_iterator
                         i(rr.products().begin()); i != rr.products().end(); ++i)
                {
                    process_ptr->registerVariableReference(
                        "_", get_variable_fullid(*i), +1);
                }
            }
            break;
        case 2:
            {
                std::stringstream ss;
                ss << "Process:/:rr" << ++num_reactions_;
                libecs::Process* const process_ptr(
                    dynamic_cast<libecs::Process*>(
                        (*model_).createEntity(
                            "DiffusionInfluencedReactionProcess",
                            libecs::FullID(ss.str()))));
                process_ptr->loadProperty("k", libecs::Polymorph(rr.k()));

                ReactionRule::reactant_container_type::const_iterator
                    r(rr.reactants().begin());
                process_ptr->registerVariableReference(
                    "_", get_variable_fullid(*r), -1);
                ++r;
                process_ptr->registerVariableReference(
                    "_", get_variable_fullid(*r), -1);
                for (ReactionRule::product_container_type::const_iterator
                         i(rr.products().begin()); i != rr.products().end(); ++i)
                {
                    process_ptr->registerVariableReference(
                        "_", get_variable_fullid(*i), +1);
                }
            }
            break;
        default:
            throw NotSupported("the number of reactants must be 1 or 2.");
        }
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

    species_container_type::const_iterator
    find_variable_fullid(const Species& sp) const
    {
        pair_first_element_unary_predicator<
            Species::serial_type, libecs::String> predicator(sp.serial());
        species_container_type::const_iterator
            i(std::find_if(
                  species_.begin(), species_.end(), predicator));
        return i;
    }

    libecs::FullID get_variable_fullid(const Species& sp) const
    {
        species_container_type::const_iterator i(find_variable_fullid(sp));
        if (i == species_.end())
        {
            throw NotFound("variable not found.");
        }
        return libecs::FullID((*i).second);
    }

    void setup_model()
    {
        (*model_).setup();
        // (*model_).setDMSearchPath("/home/kaizu/local/lib/ecell-3.2/dms");
        (*model_).setDMSearchPath("/home/kaizu/src/ecell3-spatiocyte");

        libecs::Stepper* stepper_ptr(
            (*model_).createStepper("SpatiocyteStepper", "SS"));
        stepper_ptr->setProperty(
            "VoxelRadius", libecs::Polymorph(voxel_radius_));

        (*model_).getSystem(libecs::SystemPath("/"))->setProperty(
            "StepperID", libecs::Polymorph("SS"));

        create_variable("Variable:/:GEOMETRY", CUBOID);

        const Position3 lengths(edge_lengths());
        create_variable("Variable:/:LENGTHX", lengths[0]);
        create_variable("Variable:/:LENGTHY", lengths[1]);
        create_variable("Variable:/:LENGTHZ", lengths[2]);

        create_variable("Variable:/:XYPLANE", PERIODIC);
        create_variable("Variable:/:XZPLANE", PERIODIC);
        create_variable("Variable:/:YZPLANE", PERIODIC);
        create_variable("Variable:/:VACANT", 0);

        (*model_).createEntity(
            "MoleculePopulateProcess", libecs::FullID("Process:/:populate"));
    }

    void create_variable(const libecs::String& fullid, const Real& value)
    {
        libecs::Entity* const entity_ptr(
            (*model_).createEntity("Variable", libecs::FullID(fullid)));
        entity_ptr->loadProperty("Value", libecs::Polymorph(value));
    }

    libecs::Process* get_molecule_populate_process() const
    {
        return dynamic_cast<libecs::Process*>(
            (*model_).getEntity(libecs::FullID("Process:/:populate")));
    }

    libecs::Variable* get_variable(const libecs::String& fullid) const
    {
        return dynamic_cast<libecs::Variable*>(
            (*model_).getEntity(libecs::FullID(fullid)));
    }

    void create_diffusion_process(const libecs::String& fullid, const Real& D)
    {
        const libecs::FullID vid(fullid);
        const libecs::FullID pid(
            libecs::EntityType("Process"), vid.getSystemPath(),
            "diffuse" + vid.getID());

        libecs::Process* const process_ptr(
            dynamic_cast<libecs::Process*>(
                (*model_).createEntity("DiffusionProcess", pid)));
        process_ptr->loadProperty("D", libecs::Polymorph(D));
        process_ptr->registerVariableReference("_", vid, 0);
    }

protected:

    libecs::Model* model_;
    Position3 edge_lengths_;
    Real voxel_radius_;
    Real t_;
    unsigned int num_reactions_;

    species_container_type species_;
};

} // spatiocyte

} // ecell4

#endif /* __SPATIOCYTE_WORLD_HPP */
