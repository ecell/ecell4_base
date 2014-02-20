#define BOOST_TEST_MODULE "ReactionRule_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../LatticeSpace.hpp"
#include "../SerialIDGenerator.hpp"

using namespace ecell4;

BOOST_AUTO_TEST_CASE(LatticeSpace_test_constructor)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_num_species)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
    BOOST_CHECK_EQUAL(lspace.num_species(), 0);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_has_species)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
    const Species &sp = Species("TEST");
    BOOST_CHECK(!lspace.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_update_particle)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);

    SerialIDGenerator<ParticleID> sidgen;
    ParticleID id(sidgen());

    Species sp(std::string("TEST"));
    Position3 pos(2e-7, 2.7e-7, 1.5e-7);
    Real r(1.0);
    Real d(2.3);
    Particle particle(sp, pos, r, d);

    BOOST_CHECK(lspace.update_particle(id, particle));
    BOOST_CHECK(lspace.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_num_particles)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
    BOOST_CHECK_EQUAL(lspace.num_particles(), 0);

    SerialIDGenerator<ParticleID> sidgen;
    ParticleID id(sidgen());
    Species sp(std::string("TEST"));
    Position3 pos(2e-7, 2.7e-7, 1.5e-7);
    Real r(1.0);
    Real d(2.3);
    Particle particle(sp, pos, r, d);

    ParticleID a_id(sidgen());
    Species a(std::string("ANOTHER"));
    Position3 pos1(1e-7, 2e-7, 0);
    Real r1(1.1);
    Real d1(4.3);
    Particle another(a, pos1, r1, d1);

    BOOST_CHECK(lspace.update_particle(id, particle));
    BOOST_CHECK(lspace.update_particle(a_id, another));
    BOOST_CHECK_EQUAL(lspace.num_particles(sp), 1);
    BOOST_CHECK_EQUAL(lspace.num_particles(), 2);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_list_particles)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);

    SerialIDGenerator<ParticleID> sidgen;
    ParticleID id(sidgen());
    Species sp(std::string("TEST"));
    Position3 pos(2e-7, 2.7e-7, 1.5e-7);
    Real r(1.0);
    Real d(2.3);
    Particle particle(sp, pos, r, d);

    ParticleID a_id(sidgen());
    Species a(std::string("ANOTHER"));
    Position3 pos1(1e-7, 2e-7, 0);
    Real r1(1.1);
    Real d1(4.3);
    Particle another(a, pos1, r1, d1);

    BOOST_CHECK(lspace.update_particle(id, particle));
    BOOST_CHECK(lspace.update_particle(a_id, another));

    typedef std::vector<std::pair<ParticleID, Particle> > vector;

    vector test_list(lspace.list_particles(sp));
    vector list(lspace.list_particles());
    BOOST_CHECK_EQUAL(list.size(), 2);
    BOOST_CHECK_EQUAL(test_list.size(), 1);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_add_species)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);

    Species sp(std::string("TEST"));

    BOOST_CHECK(lspace.add_species(sp));
    BOOST_CHECK(lspace.has_species(sp));

    std::vector<Species> list;
    list.push_back(sp);

    BOOST_CHECK(list == lspace.list_species());
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_add_molecule)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
    SerialIDGenerator<ParticleID> sidgen;

    Species sp(std::string("TEST"));
    BOOST_CHECK(lspace.add_species(sp));

    Coord coord(lspace.global2coord(Global(3,4,5)));
    ParticleID pid(sidgen());
    BOOST_CHECK(lspace.add_molecule(sp, coord, pid));
    BOOST_CHECK_EQUAL(lspace.num_particles(sp), 1);

    MolecularTypeBase* mt(lspace.get_molecular_type(coord));
    BOOST_CHECK(!mt->is_vacant());
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_move)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
    SerialIDGenerator<ParticleID> sidgen;

    Species sp(std::string("TEST"));
    BOOST_CHECK(lspace.add_species(sp));

    Coord coord(lspace.global2coord(Global(3,4,5)));
    ParticleID pid(sidgen());
    BOOST_CHECK(lspace.add_molecule(sp, coord, pid));

    Coord to_coord(lspace.global2coord(Global(3,5,5)));
    BOOST_CHECK(lspace.move(coord, to_coord));

    MolecularTypeBase* mt(lspace.get_molecular_type(to_coord));
    BOOST_CHECK(!mt->is_vacant());

    BOOST_CHECK(!lspace.move(coord, to_coord));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_save1)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
    SerialIDGenerator<ParticleID> sidgen;

    std::string D("1e-12"), radius("2.5e-9");
    Species sp1(std::string("A"), radius, D);
    BOOST_CHECK(lspace.add_species(sp1));

    for (int i(0); i < lspace.num_col(); i += 40)
        for (int j(0); j < lspace.num_row(); j += 40)
            for (int k(0); k < lspace.num_layer(); k += 40)
            {
                Coord coord(lspace.global2coord(Global(i, j, k)));
                ParticleID pid(sidgen());
                BOOST_CHECK(lspace.add_molecule(sp1, coord, pid));
            }

    H5::H5File fout("data.h5", H5F_ACC_TRUNC);
    const std::string hdf5path("/");
    lspace.save(&fout, hdf5path);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_save2)
{
    Position3 edge_lengths(1e-6,1e-6,1e-6);
    LatticeSpace lspace(edge_lengths);
    SerialIDGenerator<ParticleID> sidgen;

    std::string D("1e-12"), radius("2.5e-9");
    Species sp1(std::string("A"), radius, D);
    BOOST_CHECK(lspace.add_species(sp1));

    for (int i(0); i < lspace.num_col() * lspace.num_row() * lspace.num_layer(); i += 3002)
    {
        ParticleID pid(sidgen());
        BOOST_CHECK(lspace.add_molecule(sp1, i, pid));
    }

    H5::H5File fout("data_hcp.h5", H5F_ACC_TRUNC);
    const std::string hdf5path("/");
    lspace.save(&fout, hdf5path);
}
