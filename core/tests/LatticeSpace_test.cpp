#define BOOST_TEST_MODULE "ReactionRule_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../LatticeSpace.hpp"
#include "../SerialIDGenerator.hpp"

using namespace ecell4;

struct Fixture
{
    const Position3 edge_lengths;
    const Real voxel_radius;
    LatticeSpace space;
    SerialIDGenerator<ParticleID> sidgen;
    const std::string D, radius;
    const Species sp;
    Fixture() :
        edge_lengths(1e-6, 1e-6, 1e-6),
        voxel_radius(2.5e-9),
        space(edge_lengths, voxel_radius, false),
        sidgen(), D("1e-12"), radius("2.5e-9"),
        sp("A", radius, D)
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(suite, Fixture)

BOOST_AUTO_TEST_CASE(LatticeSpace_test_constructor)
{
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_num_species)
{
    BOOST_CHECK_EQUAL(space.num_species(), 0);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_has_species)
{
    BOOST_CHECK(!space.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_update_particle)
{
    ParticleID id(sidgen());

    Position3 pos(2e-7, 2.7e-7, 1.5e-7);
    Real r(1.0);
    Real d(2.3);
    Particle particle(sp, pos, r, d);

    BOOST_CHECK(space.update_particle(id, particle));
    BOOST_CHECK(space.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_num_particles)
{
    ParticleID id(sidgen());
    Position3 pos(2e-7, 2.7e-7, 1.5e-7);
    Real r(1.0), d(2.3);
    Particle particle(sp, pos, r, d);

    ParticleID a_id(sidgen());
    Species a(std::string("ANOTHER"));
    Position3 pos1(1e-7, 2e-7, 0);
    Real r1(1.1);
    Real d1(4.3);
    Particle another(a, pos1, r1, d1);

    BOOST_CHECK(space.update_particle(id, particle));
    BOOST_CHECK(space.update_particle(a_id, another));
    BOOST_CHECK_EQUAL(space.num_particles(sp), 1);
    BOOST_CHECK_EQUAL(space.num_particles(), 2);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_list_particles)
{
    ParticleID id(sidgen());
    Position3 pos(2e-7, 2.7e-7, 1.5e-7);
    Real r(1.0), d(2.3);
    Particle particle(sp, pos, r, d);

    ParticleID a_id(sidgen());
    Species a(std::string("ANOTHER"));
    Position3 pos1(1e-7, 2e-7, 0);
    Real r1(1.1);
    Real d1(4.3);
    Particle another(a, pos1, r1, d1);

    BOOST_CHECK(space.update_particle(id, particle));
    BOOST_CHECK(space.update_particle(a_id, another));

    typedef std::vector<std::pair<ParticleID, Particle> > vector;

    vector test_list(space.list_particles(sp));
    vector list(space.list_particles());
    BOOST_CHECK_EQUAL(list.size(), 2);
    BOOST_CHECK_EQUAL(test_list.size(), 1);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_add_species)
{
    BOOST_CHECK(space.add_species(sp));
    BOOST_CHECK(space.has_species(sp));

    std::vector<Species> list;
    list.push_back(sp);

    BOOST_CHECK(list == space.list_species());
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_add_remove_molecule)
{
    BOOST_CHECK(space.add_species(sp));

    Coord coord(space.global2coord(Global(3,4,5)));
    ParticleID pid(sidgen());
    BOOST_CHECK(space.add_molecule(sp, coord, pid));
    BOOST_CHECK_EQUAL(space.num_particles(sp), 1);

    MolecularTypeBase* mt(space.get_molecular_type(coord));
    BOOST_CHECK(!mt->is_vacant());

    BOOST_CHECK(space.remove_molecule(coord));
    MolecularTypeBase* vacant(space.get_molecular_type(coord));
    BOOST_CHECK(vacant->is_vacant());
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_move)
{
    BOOST_CHECK(space.add_species(sp));

    Coord coord(space.global2coord(Global(3,4,5)));
    ParticleID pid(sidgen());
    BOOST_CHECK(space.add_molecule(sp, coord, pid));

    Coord to_coord(space.global2coord(Global(3,5,5)));
    BOOST_CHECK(space.move(coord, to_coord).second);

    MolecularTypeBase* mt(space.get_molecular_type(to_coord));
    BOOST_CHECK(!mt->is_vacant());

    BOOST_CHECK(space.add_molecule(sp, coord, sidgen()));
    BOOST_CHECK(!space.move(coord, to_coord).second);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_update_molecule)
{
    Species reactant(std::string("Reactant")),
            product(std::string("Product"));
    BOOST_CHECK(space.add_species(reactant));

    Coord coord(space.global2coord(Global(3,4,5)));
    ParticleID pid(sidgen());
    BOOST_CHECK(space.add_molecule(reactant, coord, pid));

    BOOST_CHECK(space.update_molecule(coord, product));

    MolecularTypeBase* mt(space.get_molecular_type(coord));
    BOOST_ASSERT(mt->species() == product);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_save1)
{
    BOOST_CHECK(space.add_species(sp));

    for (int i(0); i < space.col_size(); i += 40)
        for (int j(0); j < space.row_size(); j += 40)
            for (int k(0); k < space.layer_size(); k += 40)
            {
                Coord coord(space.global2coord(Global(i, j, k)));
                ParticleID pid(sidgen());
                BOOST_CHECK(space.add_molecule(sp, coord, pid));
            }

    H5::H5File fout("data.h5", H5F_ACC_TRUNC);
    const std::string hdf5path("/");
    space.save(&fout, hdf5path);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_save2)
{
    BOOST_CHECK(space.add_species(sp));

    for (int i(0); i < space.col_size() * space.row_size() * space.layer_size(); i += 3002)
    {
        ParticleID pid(sidgen());
        BOOST_CHECK(space.add_molecule(sp, i, pid));
    }

    H5::H5File fout("data_hcp.h5", H5F_ACC_TRUNC);
    const std::string hdf5path("/");
    space.save(&fout, hdf5path);
}

BOOST_AUTO_TEST_SUITE_END()

struct SmallFixture
{
    const Position3 edge_lengths;
    const Real voxel_radius;
    LatticeSpace space;
    SerialIDGenerator<ParticleID> sidgen;
    const std::string D, radius;
    const Species sp;
    SmallFixture() :
        edge_lengths(2.5e-8, 2.5e-8, 2.5e-8),
        voxel_radius(2.5e-9),
        space(edge_lengths, voxel_radius, false),
        sidgen(), D("1e-12"), radius("2.5e-9"),
        sp(std::string("A"), radius, D)
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(small_suite, SmallFixture)

BOOST_AUTO_TEST_CASE(LatticeSpace_test_lattice_structure)
{
    BOOST_CHECK(space.add_species(sp));

    for (int i(0); i < space.size(); ++i)
    {
        ParticleID pid(sidgen());
        BOOST_CHECK(space.add_molecule(sp, i, pid));
    }

    H5::H5File fout("data_structure.h5", H5F_ACC_TRUNC);
    const std::string hdf5path("/");
    space.save(&fout, hdf5path);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_neighbor)
{
    BOOST_CHECK(space.add_species(sp));
    ParticleID pid(sidgen());
    BOOST_CHECK(space.add_molecule(sp, 0, pid));

    std::cout << "<<col: " << space.col_size() << ", row: "
        << space.row_size() << ", layer: " << space.layer_size() << ">>";
    for (int i(0); i < space.size(); ++i)
    {
        for (int j(0); j < 12; ++j)
        {
            const std::vector<Coord> coords(space.list_coords(sp));
            BOOST_ASSERT(coords.size() == 1);
            const Coord old(coords.at(0));
            space.move(old, i);
            const std::vector<std::pair<ParticleID, Particle> > particles(
                space.list_particles(sp));
            BOOST_ASSERT(particles.size() == 1);
            const Particle origin(particles.at(0).second);
            space.move_to_neighbor(i, j);
            const Particle neighbor(particles.at(0).second);
            const Real d(length(origin.position()-neighbor.position()));
            if (d >= 5.1e-9)
                std::cout << " [i: " << i << ", j: " << j << ", d: " << d << "] ";
            BOOST_ASSERT(d <= 5.1e-9);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

struct PeriodicFixture
{
    const Position3 edge_lengths;
    const Real voxel_radius;
    LatticeSpace space;
    SerialIDGenerator<ParticleID> sidgen;
    const std::string D, radius;
    const Species sp;
    PeriodicFixture() :
        edge_lengths(2.5e-8, 2.5e-8, 2.5e-8),
        voxel_radius(2.5e-9),
        space(edge_lengths, voxel_radius, true),
        sidgen(), D("1e-12"), radius("2.5e-9"),
        sp(std::string("A"), radius, D)
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(periodic_suite, PeriodicFixture)

BOOST_AUTO_TEST_CASE(LatticeSpace_test_periodic_col)
{
    const int col_size(space.col_size()),
              row_size(space.row_size()),
              layer_size(space.layer_size());
    const std::string hdf5path("/");
    H5::H5File fout;
    for (int i(0); i < row_size; ++i)
        for (int j(0); j < layer_size; ++j)
        {
            const Coord coord(space.global2coord(Global(0, i, j)));
            BOOST_CHECK(space.add_molecule(sp, coord, sidgen()));
        }
    fout = H5::H5File("periodic_col_0.h5", H5F_ACC_TRUNC);
    space.save(&fout, hdf5path);
    // from 0 to col_size-1
    for (int i(0); i < row_size; ++i)
        for (int j(0); j < layer_size; ++j)
        {
            const Coord coord(space.global2coord(Global(0, i, j)));
            const Integer nrnd((j&1)==1?2:3);
            //const Integer nrnd(2);
            std::pair<Coord, bool> retval(space.move_to_neighbor(coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.coord2global(retval.first).col, col_size-1);
        }
    fout = H5::H5File("periodic_col_1.h5", H5F_ACC_TRUNC);
    space.save(&fout, hdf5path);
    // from col_size-1 to 0
    for (int i(0); i < row_size; ++i)
        for (int j(0); j < layer_size; ++j)
        {
            const Coord coord(space.global2coord(Global(col_size-1, i, j)));
            const Integer nrnd((j&1)==1?4:5);
            //const Integer nrnd(4);
            std::pair<Coord, bool> retval(space.move_to_neighbor(coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.coord2global(retval.first).col, 0);
        }
    fout = H5::H5File("periodic_col_2.h5", H5F_ACC_TRUNC);
    space.save(&fout, hdf5path);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_periodic_row)
{
    const int col_size(space.col_size()),
              row_size(space.row_size()),
              layer_size(space.layer_size());
    const std::string hdf5path("/");
    H5::H5File fout;
    for (int layer(0); layer < layer_size; ++layer)
        for (int col(0); col < col_size; ++col)
        {
            const Coord coord(space.global2coord(Global(col, 0, layer)));
            BOOST_CHECK(space.add_molecule(sp, coord, sidgen()));
        }
    fout = H5::H5File("periodic_row_0.h5", H5F_ACC_TRUNC);
    space.save(&fout, hdf5path);
    // from 0 to row_size-1
    int row(0);
    for (int layer(0); layer < layer_size; ++layer)
        for (int col(0); col < col_size; ++col)
        {
            const Coord coord(space.global2coord(Global(col, row, layer)));
            const Integer nrnd(0);
            std::pair<Coord, bool> retval(space.move_to_neighbor(coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.coord2global(retval.first).row, row_size-1);
        }
    fout = H5::H5File("periodic_row_1.h5", H5F_ACC_TRUNC);
    space.save(&fout, hdf5path);
    // from row_size-1 to 0
    row = row_size - 1;
    for (int layer(0); layer < layer_size; ++layer)
        for (int col(0); col < col_size; ++col)
        {
            const Coord coord(space.global2coord(Global(col, row, layer)));
            const Integer nrnd(1);
            std::pair<Coord, bool> retval(space.move_to_neighbor(coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.coord2global(retval.first).row, 0);
        }
    fout = H5::H5File("periodic_row_2.h5", H5F_ACC_TRUNC);
    space.save(&fout, hdf5path);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_periodic_layer)
{
    const int col_size(space.col_size()),
              row_size(space.row_size()),
              layer_size(space.layer_size());
    const std::string hdf5path("/");
    H5::H5File fout;
    int layer(0);
    for (int row(0); row < row_size; ++row)
        for (int col(0); col < col_size; ++col)
        {
            const Coord coord(space.global2coord(Global(col, row, layer)));
            BOOST_CHECK(space.add_molecule(sp, coord, sidgen()));
        }
    fout = H5::H5File("periodic_layer_0.h5", H5F_ACC_TRUNC);
    space.save(&fout, hdf5path);
    // from 0 to layer_size-1
    for (int row(0); row < row_size; ++row)
        for (int col(0); col < col_size; ++col)
        {
            const Coord coord(space.global2coord(Global(col, row, layer)));
            const Integer nrnd((col&1)==1?6:7);
            std::pair<Coord, bool> retval(space.move_to_neighbor(coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.coord2global(retval.first).layer, layer_size-1);
        }
    fout = H5::H5File("periodic_layer_1.h5", H5F_ACC_TRUNC);
    space.save(&fout, hdf5path);
    return;
    // from layer_size-1 to 0
    layer = layer_size - 1;
    for (int row(0); row < row_size; ++row)
        for (int col(0); col < col_size; ++col)
        {
            const Coord coord(space.global2coord(Global(col, row, layer)));
            const Integer nrnd((col&1)==1?6:7);
            std::pair<Coord, bool> retval(space.move_to_neighbor(coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.coord2global(retval.first).layer, 0);
        }
    fout = H5::H5File("periodic_layer_2.h5", H5F_ACC_TRUNC);
    space.save(&fout, hdf5path);
}

BOOST_AUTO_TEST_SUITE_END()
