#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "BDSimulator"

#include <boost/test/included/unit_test.hpp>

#include "BDSimulator.hpp"

template<typename Tworld_, typename Trng_, typename Tpid_list_>
void inject_particles(Tworld_& world, Trng_& rng, Tpid_list_& pid_list, typename Tworld_::species_id_type const& sid, int n)
{
    typedef typename Tworld_::particle_id_pair_list particle_id_pair_list;
    typedef typename Tworld_::particle_id_pair particle_id_pair;
    typedef typename Tworld_::length_type length_type;
    typedef typename Tworld_::sphere_type sphere_type;
    typedef typename Tworld_::position_type position_type;
    typedef typename Tworld_::species_type species_type;
    boost::uniform_real<Real> ur(0., world.world_size());

    species_type const& s(world.get_species(sid));
 
    for (int i = 0; i < n; ++i)
    {
        sphere_type p(position_type(), s.radius());

        for (;;)
        {
            p.position() = position_type(ur(rng), ur(rng), ur(rng));
            if (boost::scoped_ptr<particle_id_pair_list>(
                world.check_overlap(p)) == 0)
            {
                break;
            }
            std::cerr << i << "th particle rejected" << std::endl;
        }
        
        particle_id_pair i(world.new_particle(sid, p.position()));
        pid_list.push_back(i.first);
    }
}

BOOST_AUTO_TEST_CASE(instantiation1)
{
    Model m;
    BDSimulator<BDSimulatorTraits> bds(m);
    bds.initialize();
}

BOOST_AUTO_TEST_CASE(initialization2)
{
    Model m;
    m["world_size"] = "2.0";
    m["matrix_size"] = "13";
    Model::species_type_type& S0(m.new_species_type());
    S0["name"] = "S0";
    S0["D"] = "1e-11";
    S0["radius"] = "5e-8";
    Model::species_type_type& S1(m.new_species_type());
    S1["name"] = "S1";
    S1["D"] = "0";
    S1["radius"] = "1e-8";
    Model::species_type_type& S2(m.new_species_type());
    S2["name"] = "S2";
    S2["D"] = "2e-11";
    S2["radius"] = "5e-8";
    BDSimulator<BDSimulatorTraits> bds(m);
    bds.initialize();
    BOOST_CHECK_EQUAL(2.0, bds.get_world().world_size());
    BOOST_CHECK_EQUAL(13, bds.get_world().matrix_size());
    BOOST_CHECK_EQUAL(3, bds.get_world().get_species().size());
    {
        BDSimulatorTraits::world_type::species_type const& s(
                bds.get_world().get_species(S0.id()));
        BOOST_CHECK_EQUAL(1e-11, s.D());
        BOOST_CHECK_EQUAL(5e-8, s.radius());
    }
    {
        BDSimulatorTraits::world_type::species_type const& s(
                bds.get_world().get_species(S1.id()));
        BOOST_CHECK_EQUAL(0, s.D());
        BOOST_CHECK_EQUAL(1e-8, s.radius());
    }
    {
        BDSimulatorTraits::world_type::species_type const& s(
                bds.get_world().get_species(S2.id()));
        BOOST_CHECK_EQUAL(2e-11, s.D());
        BOOST_CHECK_EQUAL(5e-8, s.radius());
    }
}
