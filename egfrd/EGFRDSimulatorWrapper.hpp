#ifndef __ECELL4_EGFRD_EGFRD_SIMULATOR_WRAPPER_HPP
#define __ECELL4_EGFRD_EGFRD_SIMULATOR_WRAPPER_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

// epdp
#include "SpeciesTypeID.hpp"
#include "ReactionRule.hpp"
#include "EGFRDSimulator.hpp"
// epdp

#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Simulator.hpp>

#include "EGFRDWorld.hpp"

#include <H5Cpp.h>
#include <hdf5.h>


namespace ecell4
{

namespace egfrd
{

class EGFRDSimulatorWrapper
    : public Simulator
{
public:

    typedef EGFRDWorld::simulator_type simulator_type;

public:

    EGFRDSimulatorWrapper(
        boost::shared_ptr<NetworkModel> model,
        boost::shared_ptr<EGFRDWorld> world,
        Integer dissociation_retry_moves = 3)
        : model_(model), world_(world)
    {
        // set the log level for epdp as L_WARNING.
        ::LoggerManager::register_logger_manager(
            "ecell.EGFRDSimulator",
            boost::shared_ptr< ::LoggerManager>(
                new ::LoggerManager("dummy", ::Logger::L_WARNING)));

        const NetworkModel::species_container_type&
            species((*model_).species());
        for (NetworkModel::species_container_type::const_iterator
                 i(species.begin()); i != species.end(); ++i)
        {
            if (!(*world_).has_species(*i))
            {
                (*world_).add_species(*i);
            }
        }

        for (NetworkModel::species_container_type::const_iterator
                 i(species.begin()); i != species.end(); ++i)
        {
            for (NetworkModel::species_container_type::const_iterator
                     j(i); j != species.end(); ++j)
            {
                if ((*model_).query_reaction_rules(*i, *j).size() == 0)
                {
                    (*world_).add_reaction_rule(
                        create_repulsive_reaction_rule(*i, *j));
                }
            }
        }

        const NetworkModel::reaction_rule_container_type&
            reaction_rules((*model_).reaction_rules());
        for (NetworkModel::reaction_rule_container_type::const_iterator
                 i(reaction_rules.begin()); i != reaction_rules.end(); ++i)
        {
            (*world_).add_reaction_rule(*i);
        }

        sim_ = boost::shared_ptr<simulator_type>(
            (*world_).create_simulator(dissociation_retry_moves));

        initialize();
    }

    // SimulatorTraits

    Real t() const
    {
        return (*sim_).t();
    }

    Real dt() const
    {
        return (*sim_).dt();
    }

    Integer num_steps() const
    {
        return (*sim_).num_steps();
    }

    void step();
    bool step(const Real& upto);

    // Optional members

    void set_t(const Real& t)
    {
        throw NotImplemented("Not implemented yet.");
    }

    void initialize()
    {
        (*sim_).initialize();
    }

    boost::shared_ptr<simulator_type> simulator() const
    {
        return sim_;
    }

	void save_hdf5_init(std::string);
	void save_hdf5(void);

protected:

    boost::shared_ptr<NetworkModel> model_;
    boost::shared_ptr<EGFRDWorld> world_;

    boost::shared_ptr<simulator_type> sim_;


	H5::H5File *file_;
	typedef struct h5_particles {
		int h5_particle_id_lot;
		int h5_particle_id_serial;
		double h5_particle_position[3];
	} h5_particles;

	typedef struct h5_particles_index {
		int h5_particle_id_serial;
		char h5_particle_name[32];

		double h5_particle_radius;
		double h5_particle_D;
	} h5_particles_index;
};

} // egfrd

} // ecell4

#endif /* __ECELL4_EGFRD_EGFRD_SIMULATOR_WRAPPER_HPP */
