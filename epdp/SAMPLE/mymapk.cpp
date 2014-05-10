// vim: foldmethod=marker
// Sakamoto EPDP Sample

#include <ecell4/egfrd_impl/config.h>
#include <stdexcept>
#include <vector>
#include <string>
#include <numeric>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <cstdlib>
#include <gsl/gsl_roots.h>

#include <boost/format.hpp>

// epdp headers
#include <ecell4/egfrd_impl/utils/range.hpp>
#include <ecell4/egfrd_impl/World.hpp>
//#include <ecell4/egfrd_impl/ParticleModel.hpp>
#include <ecell4/egfrd_impl/SpeciesType.hpp>
//#include <ecell4/egfrd_impl/SpeciesTypeID.hpp>
#include <ecell4/egfrd_impl/CuboidalRegion.hpp>
//#include <ecell4/egfrd_impl/NetworkRules.hpp>
//#include <ecell4/egfrd_impl/ReactionRule.hpp>
#include <ecell4/egfrd_impl/EGFRDSimulator.hpp>
#include <ecell4/egfrd_impl/NetworkRulesAdapter.hpp>
//#include <ecell4/egfrd_impl/GSLRandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/Model.hpp>

typedef double Real;

// Class to memize the positions of each particles
template <typename TPos_type>
class TemporaryParticleContainer {
// {{{
public:
    typedef Real radius_type;
    typedef TPos_type position_type;
    typedef std::vector< std::pair<radius_type, position_type> > particle_position_container_type;

public:
    void add(radius_type r, position_type pos)
    {
        this->particle_position_container_.push_back(
                typename particle_position_container_type::value_type(r, pos) );
    }

    particle_position_container_type
    list_particles_within_radius(radius_type r, position_type pos)
    {
        particle_position_container_type ret;
        for(    typename particle_position_container_type::iterator it = particle_position_container_.begin(); 
                it != this->particle_position_container_.end(); 
                it++) 
        {
            double radius_new(r);
            double radius_st ( it->first );
            if (this->distance(it->second, pos) < (radius_st + radius_new) ) {
                ret.push_back( *it );
            }
        }
        return ret;
    }
    double distance(position_type p1, position_type p2)
    {
        double dsq = 0.0;
        position_type sq(gsl_pow_2(p2[0] - p1[0]), gsl_pow_2(p2[1] - p1[1]), gsl_pow_2(p2[2] - p2[2]));
        return sqrt( std::accumulate(sq.begin(), sq.end(), 0.0) );
    }

private:
    particle_position_container_type particle_position_container_;
};
// }}}


int main(int argc, char **argv)
{
    // Traits typedefs  
    // {{{
    typedef ::World< ::CyclicWorldTraits<Real> > world_type;
    typedef EGFRDSimulator< ::EGFRDSimulatorTraitsBase<world_type> > simulator_type;
    typedef simulator_type::traits_type::network_rules_type network_rules_type;
    typedef simulator_type::multi_type multi_type;

    typedef ::CuboidalRegion<simulator_type::traits_type> cuboidal_region_type;
    typedef world_type::traits_type::structure_id_type structure_id_type;
    // }}}

    // Constants    
    // {{{
    const Real world_size(1e-6);
    const Integer matrix_size(3);
    const Real volume( world_size * world_size * world_size);
    const Integer N(60);
    const Real kd(0.1), U(0.5);
    const Real ka(kd * volume * (1 - U) / (U * U * N));
    const Real k2(ka), k1(kd);
    const Integer dissociation_retry_moves(3);
    // }}}

    boost::shared_ptr<ecell4::NetworkModel> ecell4_nw_model(new ecell4::NetworkModel());
    // World Definition
    // {{{
    boost::shared_ptr<world_type> world(new world_type(world_size, matrix_size));
    world_type::position_type edge_length(world_size, world_size, world_size);
    world_type::position_type pos(world_size / 2, world_size / 2, world_size / 2);

    boost::shared_ptr<cuboidal_region_type> cuboidal_region
        (new cuboidal_region_type("world", cuboidal_region_type::shape_type(pos, pos)));

    world->add_structure(cuboidal_region );
    // }}}

    // Random Number Generator (Instanciate and Initialize)
    // {{{
    boost::shared_ptr<ecell4::GSLRandomNumberGenerator> rng(new ecell4::GSLRandomNumberGenerator());
    //rng->seed(time(NULL) );
    rng->seed((unsigned long int) 0);
    world_type::traits_type::rng_type internal_rng = world_type::traits_type::rng_type( rng->handle() );
    // }}}


    // add ::SpeciesType to ::ParticleModel 
    // {{{
    ecell4::Species sp1(std::string("A"), std::string("2.5e-09"), std::string("1e-12"));
    ecell4_nw_model->add_species_attribute(sp1);

    ecell4::Species sp2(std::string("B"), std::string("2.5e-09"), std::string("1e-12"));
    ecell4_nw_model->add_species_attribute(sp2);

    ecell4::Species sp3(std::string("C"), std::string("2.5e-09"), std::string("1e-12"));
    ecell4_nw_model->add_species_attribute(sp3);

    // }}}

    // ReactionRules    
    // {{{
    // A -> B + C   k1
    // {{{
    ecell4::ReactionRule rr1( ecell4::create_unbinding_reaction_rule(sp1, sp2, sp3, k1) );
    ecell4_nw_model->add_reaction_rule(rr1);
    // }}}

    // B + C -> A   k2
    // {{{
    ecell4::ReactionRule rr2( ecell4::create_binding_reaction_rule(sp2, sp3, sp1, k2) );
    ecell4_nw_model->add_reaction_rule(rr2);
    // }}}
    // }}}

    // add ::SpeciesInfo to ::World 
    // {{{
    world->add_species(sp1.serial(), world->get_molecule_info(sp1));
    world->add_species(sp2.serial(), world->get_molecule_info(sp2));
    world->add_species(sp3.serial(), world->get_molecule_info(sp3));
    // }}}

    // Thorow particles into world at random 
    // {{{
    int number_of_particles_A(N);
    TemporaryParticleContainer<world_type::position_type> container;
    for (int cnt = 0; cnt < number_of_particles_A; cnt++) {
        // add particles at random.
        for(;;) {
            world_type::position_type particle_pos( 
                    rng->uniform(0.0, edge_length[0]), 
                    rng->uniform(0.0, edge_length[1]), 
                    rng->uniform(0.0, edge_length[2]) );
            double radius(boost::lexical_cast<double>( sp1.get_attribute("radius") ));
            if (container.list_particles_within_radius(radius, particle_pos).size() == 0) {
                std::cout << "(" << particle_pos[0] << particle_pos[1] << particle_pos[2] << ")" << std::endl;
                container.add(radius, particle_pos);
                world->new_particle( sp1.name() , particle_pos);
                break;
            }
        }
    }
    // }}}

    // world::set_all_repusive() equality section   
    // {{{
    BOOST_FOREACH( ecell4::Species temp_sp1, ecell4_nw_model->list_species() ) {
        BOOST_FOREACH( ecell4::Species temp_sp2, ecell4_nw_model->list_species() ) {
            std::vector<ecell4::ReactionRule> rrv(ecell4_nw_model->query_reaction_rules(temp_sp1, temp_sp2));
            if (rrv.size() == 0)
            {
                ecell4::ReactionRule new_reaction;
                new_reaction.add_reactant(temp_sp1);
                new_reaction.add_reactant(temp_sp2);
                new_reaction.set_k( double(0.0) );
                ecell4_nw_model->add_reaction_rule( new_reaction );
            }
        }
    }
    // }}}

    // Logger Settings 
    // {{{
    boost::shared_ptr< ::LoggerManager> logger_mng(new ::LoggerManager("dummy", ::Logger::L_WARNING));
    ::LoggerManager::register_logger_manager(
            "ecell.EGFRDSimulator",
            logger_mng
            );
    // }}}

    // EGFRDSimulator instance generated 
    // {{{
    boost::shared_ptr<network_rules_type> nw_rules_adapter(
            new network_rules_type(ecell4_nw_model));
    boost::shared_ptr< simulator_type> sim( 
            new simulator_type(
                world, 
                nw_rules_adapter,
                internal_rng,
                dissociation_retry_moves
                )
            );
    sim->initialize();
    // }}}

    // Simulation Executed
    // {{{
    Integer n_st1, n_st2, n_st3;
    n_st1 = world->get_particle_ids(sp1.name()).size();
    n_st2 = world->get_particle_ids(sp2.name()).size();
    n_st3 = world->get_particle_ids(sp3.name()).size();
    std::cout << sim->t() << "\t"
        << n_st1 << "\t"
        << n_st2 << "\t"
        << n_st3 << "\t"
        << std::endl;
    Real next_time(0.0), dt(0.02);
    for(int i(0); i < 10; i++) {
        next_time += dt;
        while(sim->step(next_time)){};
        n_st1 = world->get_particle_ids(sp1.name()).size();
        n_st2 = world->get_particle_ids(sp2.name()).size();
        n_st3 = world->get_particle_ids(sp3.name()).size();
        std::cout << sim->t() << "\t"
            << n_st1 << "\t"
            << n_st2 << "\t"
            << n_st3 << "\t"
            << std::endl;
    }
    // }}}

    // Statistics
    // {{{
    int num_single_steps_per_type[simulator_type::NUM_SINGLE_EVENT_KINDS];
    num_single_steps_per_type[simulator_type::SINGLE_EVENT_REACTION] = sim->num_single_steps_per_type( 
            simulator_type::SINGLE_EVENT_REACTION);
    num_single_steps_per_type[simulator_type::SINGLE_EVENT_ESCAPE] = sim->num_single_steps_per_type(
            simulator_type::SINGLE_EVENT_ESCAPE);
    std::cout << boost::format("%1%: %2% \n") % "SINGLE_EVENT_REACTION" % num_single_steps_per_type[simulator_type::SINGLE_EVENT_REACTION];
    std::cout << boost::format("%1%: %2% \n") % "SINGLE_EVENT_ESCAPE" % num_single_steps_per_type[simulator_type::SINGLE_EVENT_ESCAPE];

    std::cout << boost::format("%1%: %2% \n") % "PAIR_EVENT_SINGLE_REACTION_0" % sim->num_pair_steps_per_type(simulator_type::PAIR_EVENT_SINGLE_REACTION_0);
    std::cout << boost::format("%1%: %2% \n") % "PAIR_EVENT_SINGLE_REACTION_1" % sim->num_pair_steps_per_type(simulator_type::PAIR_EVENT_SINGLE_REACTION_1);
    std::cout << boost::format("%1%: %2% \n") % "PAIR_EVENT_COM_ESCAPE" % sim->num_pair_steps_per_type(simulator_type::PAIR_EVENT_COM_ESCAPE);
    std::cout << boost::format("%1%: %2% \n") % "PAIR_EVENT_IV_UNDETERMINED" % sim->num_pair_steps_per_type(simulator_type::PAIR_EVENT_IV_UNDETERMINED);
    std::cout << boost::format("%1%: %2% \n") % "PAIR_EVENT_IV_ESCAPE" % sim->num_pair_steps_per_type(simulator_type::PAIR_EVENT_IV_ESCAPE);
    std::cout << boost::format("%1%: %2% \n") % "PAIR_EVENT_IV_REACTION" % sim->num_pair_steps_per_type(simulator_type::PAIR_EVENT_IV_REACTION);

    std::cout << boost::format("%1%: %2% \n") % "NONE" % sim->num_multi_steps_per_type(multi_type::NONE);
    std::cout << boost::format("%1%: %2% \n") % "ESCAPE" % sim->num_multi_steps_per_type(multi_type::ESCAPE);
    std::cout << boost::format("%1%: %2% \n") % "REACTION" % sim->num_multi_steps_per_type(multi_type::REACTION);
    // }}}
    return 0;
}
