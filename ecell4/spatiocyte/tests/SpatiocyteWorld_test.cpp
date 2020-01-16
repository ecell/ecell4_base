#define BOOST_TEST_MODULE "SpatiocyteWorld_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include "../SpatiocyteWorld.hpp"
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Sphere.hpp>
#include <fstream>

using namespace ecell4;
using namespace ecell4::spatiocyte;

struct Fixture
{
    const Real3 edge_lengths;
    const Real voxel_radius;
    const boost::shared_ptr<GSLRandomNumberGenerator> rng;
    const boost::shared_ptr<NetworkModel> model;
    SpatiocyteWorld world;

    Fixture()
        : edge_lengths(1e-6, 1e-6, 1e-6), voxel_radius(1e-8),
          rng(new GSLRandomNumberGenerator()), model(new NetworkModel),
          world(edge_lengths, voxel_radius, rng)
    {
        world.bind_to(model);
    }
};

BOOST_FIXTURE_TEST_SUITE(suite, Fixture)

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_constructor) {}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_t)
{
    BOOST_CHECK_EQUAL(world.t(), 0);
    world.set_t(23.4);
    BOOST_CHECK_EQUAL(world.t(), 23.4);
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_num_species)
{
    BOOST_CHECK_EQUAL(world.list_species().size(), 0);
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_has_species)
{
    Species sp(std::string("Species"));
    BOOST_CHECK(!world.has_species(sp));
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_list_particles)
{
    std::vector<std::pair<ParticleID, Particle>> particles(
        world.list_particles());
    BOOST_CHECK_EQUAL(particles.size(), 0);
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_update_particles)
{
    SerialIDGenerator<ParticleID> sidgen;
    ParticleID pid(sidgen());
    Species sp(std::string("A"));
    const Real3 pos(2e-7, 1e-7, 0);
    Real r(0);
    Real d(0);
    Particle p(sp, pos, r, d);

    model->add_species_attribute(sp);
    world.update_particle(pid, p);

    BOOST_CHECK(world.has_species(sp));
    BOOST_CHECK(world.has_particle(pid));
    BOOST_CHECK_EQUAL(world.list_particles().size(), 1);
    BOOST_CHECK_EQUAL(world.list_particles(sp).size(), 1);
}

// BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_register_species)
// {
//     Species sp(std::string("TEST"));
//
//     BOOST_CHECK(world.register_species(sp));
//     BOOST_CHECK(world.has_species(sp));
//
//     std::vector<Species> list;
//     list.push_back(sp);
//
//     BOOST_CHECK(list == world.list_species());
// }

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_add_molecule)
{
    const Species sp("TEST", 1e-8, 1e-12);
    model->add_species_attribute(sp);

    const Voxel voxel(world.get_voxel_nearby(edge_lengths / 2.0));
    // BOOST_CHECK(world.place_voxel(sp, coord).second);
    BOOST_CHECK(world.new_particle(sp, voxel));
    BOOST_CHECK_EQUAL(world.num_particles(sp), 1);

    boost::shared_ptr<const VoxelPool> mt(voxel.get_voxel_pool());
    BOOST_CHECK(!mt->is_vacant());
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_add_molecules)
{
    const Species sp("TEST", 1e-8, 1e-12);
    model->add_species_attribute(sp);

    const Integer N(60);
    BOOST_CHECK(world.add_molecules(sp, N));
    BOOST_CHECK_EQUAL(world.num_particles(sp), N);
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_neighbor)
{
    const Voxel voxel(world.get_voxel_nearby(edge_lengths / 2.0));
    const Real3 cp(voxel.position());

    const Species sp("TEST", 1e-8, 1e-12);
    model->add_species_attribute(sp);

    for (Integer i(0); i < world.num_neighbors(voxel); ++i)
    {
        world.new_particle(sp, world.get_neighbor(voxel, i));
    }
    std::vector<std::pair<ParticleID, Particle>> particles(
        world.list_particles());
    for (std::vector<std::pair<ParticleID, Particle>>::iterator itr(
             particles.begin());
         itr != particles.end(); ++itr)
    {
        Real3 pos((*itr).second.position());
        BOOST_ASSERT(length(pos - cp) < voxel_radius * 2.1);
    }

#ifdef WITH_HDF5
    world.save("neighbor.h5");
#endif
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_add_shape)
{
    const Species sp("TEST", 1e-8, 1e-12);
    model->add_species_attribute(sp);

    boost::shared_ptr<const Sphere> sphere(
        new Sphere(Real3(5e-7, 5e-7, 5e-7), 5e-7 * 1.5));

    const Integer n(world.add_structure(sp, sphere));
    BOOST_ASSERT(n > 0);
    BOOST_CHECK_EQUAL(world.num_particles(sp), n);

#ifdef WITH_HDF5
    world.save("sphere.h5");
#endif
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_move)
{
    const Species sp("TEST", 1e-8, 1e-12);
    model->add_species_attribute(sp);

    const Voxel from(world.get_voxel_nearby(Real3(0.3e-6, 0.5e-6, 0.5e-6)));
    const Voxel to(world.get_voxel_nearby(Real3(0.5e-6, 0.5e-6, 0.5e-6)));

    BOOST_CHECK(world.new_particle(sp, from));
    BOOST_CHECK(world.move(from, to));

    boost::shared_ptr<const VoxelPool> mt(to.get_voxel_pool());
    BOOST_CHECK(!mt->is_vacant());

    BOOST_CHECK(world.move(from, to));
}

BOOST_AUTO_TEST_CASE(SpatiocyteWorld_test_structure)
{
    const Species membrane("Membrane", 2.5e-9, 0);
    const Species sp("TEST", 1e-8, 1e-12, "Membrane");
    model->add_species_attribute(membrane);
    model->add_species_attribute(sp);

    boost::shared_ptr<const Sphere> sphere(
        new Sphere(Real3(2.5e-7, 2.5e-7, 2.5e-7), 2e-7));

    BOOST_CHECK(world.add_structure(membrane, sphere) == 5892);
    BOOST_CHECK(!world.new_particle(
        Particle(sp, Real3(2.5e-7, 2.5e-7, 4.5e-7), 2.5e-9, 1e-12)));
    BOOST_CHECK(world.new_particle(Particle(
        sp, Real3(2.5e-7, 2.5e-7, 4.5e-7 - voxel_radius * 2), 2.5e-9, 1e-12)));

#ifdef WITH_HDF5
    world.save("structure.h5");
#endif
}

BOOST_AUTO_TEST_SUITE_END()
