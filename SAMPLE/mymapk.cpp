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

// epdp headers
#include <ecell4/egfrd_impl/utils/range.hpp>
#include <ecell4/egfrd_impl/World.hpp>
#include <ecell4/egfrd_impl/ParticleModel.hpp>
#include <ecell4/egfrd_impl/SpeciesType.hpp>
#include <ecell4/egfrd_impl/SpeciesTypeID.hpp>
#include <ecell4/egfrd_impl/CuboidalRegion.hpp>
#include <ecell4/egfrd_impl/NetworkRules.hpp>
#include <ecell4/egfrd_impl/ReactionRule.hpp>
#include <ecell4/egfrd_impl/EGFRDSimulator.hpp>
#include <ecell4/egfrd_impl/NetworkRulesAdapter.hpp>
//#include <ecell4/egfrd_impl/GSLRandomNumberGenerator.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/NetworkModel.hpp>

typedef double Real;

typedef ::World< ::CyclicWorldTraits<Real, Real> > world_type;

double distance_sq(world_type::position_type p1, world_type::position_type p2)
{
    double dsq = 0.0;
    world_type::position_type sq(gsl_pow_2(p2[0] - p1[0]), gsl_pow_2(p2[1] - p1[1]), gsl_pow_2(p2[2] - p2[2]));
    return std::accumulate(sq.begin(), sq.end(), 0.0);
}

// Class to memize the positions of each particles
class TemporaryParticleContainer {  
// {{{
public:
    typedef std::vector<std::pair< boost::shared_ptr< ::SpeciesType>, world_type::position_type> >  particle_position_container;
    TemporaryParticleContainer(void) {;}

    void add( boost::shared_ptr< ::SpeciesType> st, world_type::position_type pos)
    {
        this->container_.push_back( particle_position_container::value_type(st, pos));
    }

    particle_position_container 
    list_particles_within_radius(boost::shared_ptr< ::SpeciesType> st, world_type::position_type &pos)
    {
        particle_position_container ret;
        for(particle_position_container::iterator it = container_.begin(); it != container_.end(); it++) {
            double radius_new( atof(((*st)["radius"]).c_str()) );
            double radius_st(  atof(((*(it->first))["radius"]).c_str()) );
            if (distance_sq(it->second, pos) < gsl_pow_2(radius_new) ) {
                ret.push_back( *it );
            }
        }
        return ret;
    }

private:
    particle_position_container container_;
};
// }}}  

int main(int argc, char **argv)
{
    // Traits typedefs  
    // {{{
    typedef ::World< ::CyclicWorldTraits<Real, Real> > world_type;
    typedef ::ParticleModel particle_model_type;
    typedef EGFRDSimulator< ::EGFRDSimulatorTraitsBase<world_type> > simulator_type;
    typedef simulator_type::traits_type::network_rules_type network_rules_type;

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

    // World Definition
    // {{{
    boost::shared_ptr<world_type> world(new world_type(world_size, matrix_size));
    world_type::position_type edge_length(world_size, world_size, world_size);
    world_type::position_type pos(world_size / 2, world_size / 2, world_size / 2);
    world->add_structure( boost::shared_ptr<cuboidal_region_type>(
                new cuboidal_region_type("world", cuboidal_region_type::shape_type(pos, pos))));
    // }}}

    // Random Number Generator (Instanciate and Initialize)
    // {{{
    boost::shared_ptr<ecell4::GSLRandomNumberGenerator> rng(new ecell4::GSLRandomNumberGenerator());
    particle_model_type model;
    //rng->seed( (unsigned long int)0 );
    rng->seed(time(NULL) );
    // }}}

    world_type::traits_type::rng_type internal_rng = world_type::traits_type::rng_type( rng->handle() );

    // add ::SpeciesType to ::ParticleModel 
    // {{{
    boost::shared_ptr< ::SpeciesType> st1(new ::SpeciesType());
    (*st1)["name"] = std::string("A");
    (*st1)["D"] = std::string("1e-12");
    (*st1)["radius"] = std::string("2.5e-9");
    model.add_species_type(st1);

    boost::shared_ptr< ::SpeciesType> st2(new ::SpeciesType());
    (*st2)["name"] = std::string("A");
    (*st2)["D"] = std::string("1e-12");
    (*st2)["radius"] = std::string("2.5e-9");
    model.add_species_type(st2);

    boost::shared_ptr< ::SpeciesType> st3(new ::SpeciesType());
    (*st3)["name"] = std::string("A");
    (*st3)["D"] = std::string("1e-12");
    (*st3)["radius"] = std::string("2.5e-9");
    model.add_species_type(st3);
    // }}}

    // ReactionRules    
    // {{{
    // A -> B + C   k1
    // {{{
    std::vector< ::SpeciesTypeID> products;
    products.push_back(st2->id());
    products.push_back(st3->id());
    model.network_rules().add_reaction_rule( new_reaction_rule(st1->id(), products, k1) );
    // }}}

    // B + C -> A   k2
    // {{{
    products.clear();
    products.push_back(st1->id());
    model.network_rules().add_reaction_rule( new_reaction_rule(st2->id(), st3->id(), products, k2) );
    // }}}
    // }}}

    // add ::SpeciesInfo to ::World 
    // {{{
    //  st1 {{{
    const std::string &structure_id((*st1)["structure"]);
    world->add_species( world_type::traits_type::species_type(
                st1->id(), 
                boost::lexical_cast<world_type::traits_type::D_type>( (*st1)["D"] ),
                boost::lexical_cast<world_type::length_type>( (*st1)["radius"] ),
                boost::lexical_cast<structure_id_type>( structure_id.empty() ? "world" : structure_id )));
    // }}}

    //  st2 {{{
    const std::string &structure_id2((*st2)["structure"]);
    world->add_species( world_type::traits_type::species_type(
                st2->id(), 
                boost::lexical_cast<world_type::traits_type::D_type>( (*st2)["D"] ),
                boost::lexical_cast<world_type::length_type>( (*st2)["radius"] ),
                boost::lexical_cast<structure_id_type>( structure_id.empty() ? "world" : structure_id2 )));
    // }}}

    //  st3 {{{
    const std::string &structure_id3((*st3)["structure"]);
    world->add_species( world_type::traits_type::species_type(
                st3->id(), 
                boost::lexical_cast<world_type::traits_type::D_type>( (*st3)["D"] ),
                boost::lexical_cast<world_type::length_type>( (*st3)["radius"] ),
                boost::lexical_cast<structure_id_type>( structure_id.empty() ? "world" : structure_id3 )));
    // }}}
    // }}}

    // Thorow particles into world at random 
    // {{{
    int number_of_particles_A(N);
    TemporaryParticleContainer container;
    for (int cnt = 0; cnt < number_of_particles_A; cnt++) {
        // add particles at random.
        for(;;) {
            world_type::position_type particle_pos( rng->uniform(0.0, edge_length[0]), rng->uniform(0.0, edge_length[1]), rng->uniform(0.0, edge_length[2]) );
            if (container.list_particles_within_radius(st1, particle_pos).size() == 0) {
                //std::cout << "(" << particle_pos[0] << particle_pos[1] << particle_pos[2] << ")" << std::endl;
                container.add(st1, particle_pos);
                world->new_particle(st1->id(), particle_pos);
                break;
            }
        }
    }
    // }}}

    // world::set_all_repusive() equality section   
    // {{{
    BOOST_FOREACH( boost::shared_ptr< ::SpeciesType> temp_st1, model.get_species_types()) {
        BOOST_FOREACH( boost::shared_ptr< ::SpeciesType> temp_st2, model.get_species_types()) {
            boost::scoped_ptr< ::NetworkRules::reaction_rule_generator> gen( model.network_rules().query_reaction_rule( temp_st1->id(), temp_st2->id()));
            if (!gen) {
                const::std::vector< ::SpeciesTypeID> products;
                model.network_rules().add_reaction_rule( ::new_reaction_rule(temp_st1->id()(), temp_st2->id(), products, 0.0) );
            }
        }
    }   // }}}

    // Logger Settings 
    // {{{
    ::LoggerManager::register_logger_manager(
            "ecell.EGFRDSimulator",
            boost::shared_ptr< ::LoggerManager>(
                new ::LoggerManager("dummy", ::Logger::L_WARNING)));
    // }}}

    // EGFRDSimulator instance generated 
    // {{{
    boost::shared_ptr<ecell4::NetworkModel> ecell4_nw_model(new ecell4::NetworkModel());
    boost::shared_ptr< simulator_type> sim( 
            new simulator_type(
                world, 
                //boost::shared_ptr<network_rules_type>(new network_rules_type(model.network_rules())),
                boost::shared_ptr<network_rules_type>(new network_rules_type(ecell4_nw_model) ),
                internal_rng,
                dissociation_retry_moves
                )
            );
    sim->initialize();
    // }}}

    // Simulation Executed
    // {{{
    Integer n_st1, n_st2, n_st3;
    n_st1 = world->get_particle_ids(st1->id()).size();
    n_st2 = world->get_particle_ids(st2->id()).size();
    n_st3 = world->get_particle_ids(st3->id()).size();
    std::cout << sim->t() << "\t"
        << n_st1 << "\t"
        << n_st2 << "\t"
        << n_st3 << "\t"
        << std::endl;
    Real next_time(0.0), dt(0.02);
    for(int i(0); i < 100; i++) {
        next_time += dt;
        while(sim->step(next_time)){};
        n_st1 = world->get_particle_ids(st1->id()).size();
        n_st2 = world->get_particle_ids(st2->id()).size();
        n_st3 = world->get_particle_ids(st3->id()).size();
        std::cout << sim->t() << "\t"
            << n_st1 << "\t"
            << n_st2 << "\t"
            << n_st3 << "\t"
            << std::endl;
    }
    // }}}

    return 0;
}
