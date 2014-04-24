#define BOOST_TEST_MODULE "ReactionRule_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "MolecularType.hpp"
#include "VacantType.hpp"
#include "../LatticeSpace.hpp"
#include "../SerialIDGenerator.hpp"

using namespace ecell4;

struct Fixture
{
    const Position3 edge_lengths;
    const Real voxel_radius;
    LatticeSpace space;
    SerialIDGenerator<ParticleID> sidgen;
    const Real D, radius;
    const Species sp;
    Fixture() :
        edge_lengths(2.5e-8, 2.5e-8, 2.5e-8),
        voxel_radius(2.5e-9),
        space(edge_lengths, voxel_radius, false),
        sidgen(), D(1e-12), radius(2.5e-9),
        sp("A", "2.5e-9", "1e-12")
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

    Position3 pos(2e-8, 1.7e-8, 1.5e-8);
    Real r(1.0);
    Real d(2.3);
    Particle particle(sp, pos, r, d);

    BOOST_CHECK(space.update_particle(id, particle));
    BOOST_CHECK(space.has_species(sp));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_num_particles)
{
    ParticleID id(sidgen());
    Position3 pos(2e-8, 1.7e-8, 1.5e-8);
    Real r(1.0), d(2.3);
    Particle particle(sp, pos, r, d);

    ParticleID a_id(sidgen());
    Species a(std::string("ANOTHER"));
    Position3 pos1(1e-8, 2e-8, 0);
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
    Position3 pos(2e-8, 1.7e-8, 1.5e-8);
    Real r(1.0), d(2.3);
    Particle particle(sp, pos, r, d);

    ParticleID a_id(sidgen());
    Species a(std::string("ANOTHER"));
    Position3 pos1(1e-8, 2e-8, 0);
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

// BOOST_AUTO_TEST_CASE(LatticeSpace_test_register_species)
// {
//     BOOST_CHECK(space.register_species(sp));
//     BOOST_CHECK(space.has_species(sp));
// 
//     std::vector<Species> list;
//     list.push_back(sp);
// 
//     BOOST_CHECK(list == space.list_species());
// }

/*
 * for Simulator
 */
BOOST_AUTO_TEST_CASE(LatticeSpace_test_coordinate)
{
    for (Integer col(0); col < space.col_size(); ++col)
        for (Integer row(0); row < space.row_size(); ++row)
            for (Integer layer(0); layer < space.layer_size(); ++layer)
            {
                const Global global(col, row, layer);
                const LatticeSpace::private_coordinate_type private_coord(
                        space.global2private_coord(global));
                const LatticeSpace::coordinate_type coord(space.global2coord(global));
                BOOST_CHECK_EQUAL(private_coord, space.coord2private(coord));
                BOOST_CHECK_EQUAL(space.private2coord(private_coord), coord);
            }
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_add_remove_molecule)
{
    const LatticeSpace::private_coordinate_type private_coord(
            space.global2private_coord(Global(3,4,5)));
    ParticleID pid(sidgen());
    BOOST_CHECK(space.update_voxel_private(
        pid, Voxel(sp, private_coord, radius, D)));
    BOOST_CHECK_EQUAL(space.num_particles(sp), 1);

    const MolecularTypeBase* mt(space.get_molecular_type(private_coord));
    BOOST_CHECK(!mt->is_vacant());

    BOOST_CHECK(space.remove_molecule(private_coord));
    const MolecularTypeBase* vacant(space.get_molecular_type(private_coord));
    BOOST_CHECK(vacant->is_vacant());
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_move)
{
    const Global global0(3,4,5);
    const LatticeSpace::private_coordinate_type private_coord(
            space.global2private_coord(global0));
    const LatticeSpace::coordinate_type coord(
            space.global2coord(global0));

    ParticleID pid(sidgen());
    BOOST_CHECK(space.update_voxel_private(
        pid, Voxel(sp, private_coord, radius, D)));

    MolecularTypeBase* from_mt(space.get_molecular_type(private_coord));
    BOOST_CHECK(!from_mt->is_vacant());

    const Global global1(3,5,5);
    const LatticeSpace::private_coordinate_type private_to_coord(
            space.global2private_coord(global1));
    const LatticeSpace::coordinate_type to_coord(
            space.global2coord(global1));

    BOOST_CHECK(space.move(coord, to_coord));

    MolecularTypeBase* mt(space.get_molecular_type(private_to_coord));
    BOOST_CHECK(!mt->is_vacant());

    BOOST_CHECK(space.update_voxel_private(
        sidgen(), Voxel(sp, private_coord, radius, D)));
    BOOST_CHECK(!space.move(coord, to_coord));
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_update_molecule)
{
    Species reactant(std::string("Reactant")),
            product(std::string("Product"));

    const Global global(3,4,5);
    const LatticeSpace::private_coordinate_type private_coord(
            space.global2private_coord(global));

    ParticleID pid(sidgen());
    BOOST_CHECK(space.update_voxel_private(
        pid, Voxel(reactant, private_coord, radius, D)));
    BOOST_CHECK(space.update_voxel_private(
        Voxel(product, private_coord, radius, D)));

    const MolecularTypeBase* mt(space.get_molecular_type(private_coord));
    BOOST_ASSERT(mt->species() == product);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_save)
{
    for (Integer col(0); col < space.col_size(); col += 2)
        for (Integer row(0); row < space.row_size(); row += 2)
            for (Integer layer(0); layer < space.layer_size(); layer += 2)
            {
                const LatticeSpace::private_coordinate_type private_coord(
                        space.global2private_coord(Global(col, row, layer)));
                ParticleID pid(sidgen());
                BOOST_CHECK(space.update_voxel_private(
                    pid, Voxel(sp, private_coord, radius, D)));
            }

    H5::H5File fout("data.h5", H5F_ACC_TRUNC);
    boost::scoped_ptr<H5::Group>
        group(new H5::Group(fout.createGroup("LatticeSpace")));
    space.save(group.get());
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_lattice_structure)
{
    for (LatticeSpace::coordinate_type coord(0); coord < space.size(); ++coord)
    {
        ParticleID pid(sidgen());
        const LatticeSpace::private_coordinate_type private_coord(
                space.coord2private(coord));
        BOOST_CHECK(space.update_voxel_private(
            pid, Voxel(sp, private_coord, radius, D)));
    }

    H5::H5File fout("data_structure.h5", H5F_ACC_TRUNC);
    boost::scoped_ptr<H5::Group>
        group(new H5::Group(fout.createGroup("LatticeSpace")));
    space.save(group.get());
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_neighbor)
{
    ParticleID pid(sidgen());
    BOOST_CHECK(space.update_voxel_private(
        pid, Voxel(sp, space.coord2private(0), radius, D)));

    std::cout << "<<col: " << space.col_size() << ", row: "
        << space.row_size() << ", layer: " << space.layer_size() << ">>";
    for (int coord(0); coord < space.size(); ++coord)
    {
        for (int j(0); j < 12; ++j)
        {
            const std::vector<LatticeSpace::coordinate_type> coords(
                    space.list_coords(sp));
            BOOST_ASSERT(coords.size() == 1);

            const LatticeSpace::coordinate_type old(coords.at(0));
            space.move(old, coord);

            const std::vector<std::pair<ParticleID, Particle> > particles(
                space.list_particles(sp));
            BOOST_ASSERT(particles.size() == 1);
            const Particle origin(particles.at(0).second);
            const LatticeSpace::private_coordinate_type private_coord(
                    space.coord2private(coord));
            space.move_to_neighbor(private_coord, j);
            const Particle neighbor(particles.at(0).second);
            const Real d(length(origin.position()-neighbor.position()));
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
    const Real D, radius;
    const Species sp;
    PeriodicFixture() :
        edge_lengths(2.5e-8, 2.5e-8, 2.5e-8),
        voxel_radius(2.5e-9),
        space(edge_lengths, voxel_radius, true),
        sidgen(), D(1e-12), radius(2.5e-9),
        sp(std::string("A"), "2.5e-9", "1e-12")
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(periodic_suite, PeriodicFixture)

BOOST_AUTO_TEST_CASE(LatticeSpace_test_periodic_col)
{
    std::cerr << " < periodic_col > ";
    const int col_size(space.col_size()),
              row_size(space.row_size()),
              layer_size(space.layer_size());
    H5::H5File fout;
    H5::Group group;
    for (int i(0); i < row_size; ++i)
        for (int j(0); j < layer_size; ++j)
        {
            const LatticeSpace::private_coordinate_type private_coord(
                    space.global2private_coord(Global(0, i, j)));
            BOOST_CHECK(space.update_voxel_private(
                sidgen(), Voxel(sp, private_coord, radius, D)));
        }
    fout = H5::H5File("periodic_col_0.h5", H5F_ACC_TRUNC);
    group = fout.createGroup("LatticeSpace");
    space.save(&group);

    // from 0 to col_size-1
    for (int i(0); i < row_size; ++i)
        for (int j(0); j < layer_size; ++j)
        {
            const LatticeSpace::private_coordinate_type private_coord(
                    space.global2private_coord(Global(0, i, j)));
            const Integer nrnd((j&1)==1?2:3);
            std::pair<LatticeSpace::private_coordinate_type, bool> retval(
                    space.move_to_neighbor(private_coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.private_coord2global(retval.first).col,
                    col_size-1);
        }
    fout = H5::H5File("periodic_col_1.h5", H5F_ACC_TRUNC);
    group = fout.createGroup("LatticeSpace");
    space.save(&group);

    // from col_size-1 to 0
    for (int i(0); i < row_size; ++i)
        for (int j(0); j < layer_size; ++j)
        {
            const LatticeSpace::private_coordinate_type private_coord(
                    space.global2private_coord(Global(col_size-1, i, j)));
            const Integer nrnd((j&1)==1?4:5);
            std::pair<LatticeSpace::private_coordinate_type, bool> retval(
                    space.move_to_neighbor(private_coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.private_coord2global(retval.first).col, 0);
        }
    fout = H5::H5File("periodic_col_2.h5", H5F_ACC_TRUNC);
    group = fout.createGroup("LatticeSpace");
    space.save(&group);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_periodic_row)
{
    const int col_size(space.col_size()),
              row_size(space.row_size()),
              layer_size(space.layer_size());
    H5::H5File fout;
    H5::Group group;
    for (int layer(0); layer < layer_size; ++layer)
        for (int col(0); col < col_size; ++col)
        {
            const LatticeSpace::private_coordinate_type private_coord(
                    space.global2private_coord(Global(col, 0, layer)));
            BOOST_CHECK(space.update_voxel_private(
                sidgen(), Voxel(sp, private_coord, radius, D)));
        }
    fout = H5::H5File("periodic_row_0.h5", H5F_ACC_TRUNC);
    group = fout.createGroup("LatticeSpace");
    space.save(&group);
    // from 0 to row_size-1
    int row(0);
    for (int layer(0); layer < layer_size; ++layer)
        for (int col(0); col < col_size; ++col)
        {
            const LatticeSpace::private_coordinate_type private_coord(
                    space.global2private_coord(Global(col, row, layer)));
            const Integer nrnd(0);
            std::pair<LatticeSpace::private_coordinate_type, bool> retval(
                    space.move_to_neighbor(private_coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.private_coord2global(retval.first).row,
                    row_size-1);
        }
    fout = H5::H5File("periodic_row_1.h5", H5F_ACC_TRUNC);
    group = fout.createGroup("LatticeSpace");
    space.save(&group);
    // from row_size-1 to 0
    row = row_size - 1;
    for (int layer(0); layer < layer_size; ++layer)
        for (int col(0); col < col_size; ++col)
        {
            const LatticeSpace::private_coordinate_type private_coord(
                    space.global2private_coord(Global(col, row, layer)));
            const Integer nrnd(1);
            std::pair<LatticeSpace::private_coordinate_type, bool> retval(
                    space.move_to_neighbor(private_coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.private_coord2global(retval.first).row, 0);
        }
    fout = H5::H5File("periodic_row_2.h5", H5F_ACC_TRUNC);
    group = fout.createGroup("LatticeSpace");
    space.save(&group);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_periodic_layer)
{
    const int col_size(space.col_size()),
              row_size(space.row_size()),
              layer_size(space.layer_size());
    H5::H5File fout;
    H5::Group group;
    int layer(0);
    for (int row(0); row < row_size; ++row)
        for (int col(0); col < col_size; ++col)
        {
            const LatticeSpace::private_coordinate_type private_coord(
                    space.global2private_coord(Global(col, row, layer)));
            BOOST_CHECK(space.update_voxel_private(
                sidgen(), Voxel(sp, private_coord, radius, D)));
        }
    fout = H5::H5File("periodic_layer_0.h5", H5F_ACC_TRUNC);
    group = fout.createGroup("LatticeSpace");
    space.save(&group);
    // from 0 to layer_size-1
    for (int row(0); row < row_size; ++row)
        for (int col(0); col < col_size; ++col)
        {
            const LatticeSpace::private_coordinate_type private_coord(
                    space.global2private_coord(Global(col, row, layer)));
            const Integer nrnd((col&1)==1?6:7);
            std::pair<LatticeSpace::private_coordinate_type, bool> retval(
                    space.move_to_neighbor(private_coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.private_coord2global(retval.first).layer,
                    layer_size-1);
        }
    fout = H5::H5File("periodic_layer_1.h5", H5F_ACC_TRUNC);
    group = fout.createGroup("LatticeSpace");
    space.save(&group);
    return;
    // from layer_size-1 to 0
    layer = layer_size - 1;
    for (int row(0); row < row_size; ++row)
        for (int col(0); col < col_size; ++col)
        {
            const LatticeSpace::private_coordinate_type private_coord(
                    space.global2private_coord(Global(col, row, layer)));
            const Integer nrnd((col&1)==1?6:7);
            std::pair<LatticeSpace::private_coordinate_type, bool> retval(
                    space.move_to_neighbor(private_coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.private_coord2global(retval.first).layer, 0);
        }
    fout = H5::H5File("periodic_layer_2.h5", H5F_ACC_TRUNC);
    group = fout.createGroup("LatticeSpace");
    space.save(&group);
}

BOOST_AUTO_TEST_SUITE_END()
