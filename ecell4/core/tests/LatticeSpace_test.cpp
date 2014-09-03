#define BOOST_TEST_MODULE "ReactionRule_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/MolecularType.hpp>
#include <ecell4/core/VacantType.hpp>
#include <ecell4/core/LatticeSpace.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>

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
    Position3 pos1(1e-8, 2e-8, 1e-9);
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
    Position3 pos1(1e-8, 2e-8, 1e-9);
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
    for (LatticeSpace::coordinate_type coord(0); coord < space.size(); ++coord)
    {
        const LatticeSpace::private_coordinate_type private_coord(
                space.coord2private(coord));
        BOOST_CHECK_EQUAL(coord, space.private2coord(private_coord));
    }

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

BOOST_AUTO_TEST_CASE(LatticeSpace_test_coordinate_global_translation)
{
    for (LatticeSpace::coordinate_type coord(0); coord < space.size(); ++coord)
    {
        const Global global(space.coord2global(coord));
        LatticeSpace::coordinate_type created_coord(
                space.global2coord(global));
        BOOST_CHECK_EQUAL(coord, created_coord);
    }
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_coordinate_position_translation)
{
    const Position3 private_origin(space.private2position(0));
    BOOST_ASSERT(private_origin[0] < 0);
    BOOST_ASSERT(private_origin[1] < 0);
    BOOST_ASSERT(private_origin[2] < 0);

    const LatticeSpace::private_coordinate_type origin(
                (space.col_size() + 3) * (space.row_size() + 2) + 1);
    const Position3 origin_p(space.private2position(origin));
    BOOST_ASSERT(origin_p[0] == 0);
    BOOST_ASSERT(origin_p[1] == 0);
    BOOST_ASSERT(origin_p[2] == 0);

    for (Integer i(0); i < 12; ++i)
    {
        const Position3 neighbor(
                space.private2position(space.get_neighbor(origin, i)));
        const LatticeSpace::private_coordinate_type coord(
                space.position2private(origin_p * 0.7 + neighbor * 0.3));
        BOOST_CHECK_EQUAL(origin, coord);
    }

    Integer size(
            (space.col_size()+2) * (space.layer_size() + 2) * (space.row_size() + 2));
    for (LatticeSpace::private_coordinate_type coord(0); coord < size; ++coord)
    {
        const Position3 pos(space.private2position(coord));
        const Global global(space.position2global(pos));
        const LatticeSpace::private_coordinate_type created_coord(
                space.position2private(pos));
        BOOST_CHECK_EQUAL(coord, created_coord);
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

    BOOST_CHECK(space.remove_voxel_private(private_coord));
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

BOOST_AUTO_TEST_CASE(LatticeSpace_test_update_voxel)
{
    const ParticleID pid(sidgen());
    for (LatticeSpace::coordinate_type coord(0); coord < space.size(); ++coord)
    {
        const Position3 pos(space.coordinate2position(coord));
        BOOST_CHECK(space.update_voxel(pid, Voxel(sp, coord, radius, D)));
        BOOST_CHECK_EQUAL(space.num_particles(), 1);
        std::pair<ParticleID, Particle> pair(space.list_particles()[0]);
        BOOST_CHECK_EQUAL(pid, pair.first);
        BOOST_CHECK_EQUAL(pos, pair.second.position());
        BOOST_CHECK_EQUAL(radius, pair.second.radius());
        BOOST_CHECK_EQUAL(D, pair.second.D());
        //BOOST_CHECK_EQUAL(sp, pair.second.species());
        //[TODO] Species is not comparable.
    }
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_save_and_load)
{
    const ParticleID pid(sidgen());
    for (LatticeSpace::coordinate_type coord(0); coord < space.size(); ++coord)
    {
        const Position3 pos(space.coordinate2position(coord));
        BOOST_CHECK(space.update_voxel(pid, Voxel(sp, coord, radius, D)));

        H5::H5File fout("data.h5", H5F_ACC_TRUNC);
        boost::scoped_ptr<H5::Group>
            group(new H5::Group(fout.createGroup("LatticeSpace")));
        space.save(group.get());
        fout.close();

        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File("data.h5", H5F_ACC_RDONLY));
        const H5::Group groupin(fin->openGroup("LatticeSpace"));
        LatticeSpace space2(Position3(3e-8, 3e-8, 3e-8), voxel_radius);
        space2.load(groupin);
        fin->close();

        BOOST_CHECK_EQUAL(space.edge_lengths(), space2.edge_lengths());
        std::pair<ParticleID, Particle> pair(space2.list_particles()[0]);
        BOOST_CHECK_EQUAL(pid, pair.first);
        BOOST_CHECK_EQUAL(pos, pair.second.position());
        BOOST_CHECK_EQUAL(radius, pair.second.radius());
        BOOST_CHECK_EQUAL(D, pair.second.D());
    }
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
    for (LatticeSpace::coordinate_type coord(0); coord < space.size(); ++coord)
    {
        Position3 center(space.coordinate2position(coord));
        for (int i(0); i < 12; ++i)
        {
            LatticeSpace::private_coordinate_type neighbor(
                    space.get_neighbor(space.coord2private(coord), i));
            if (!space.is_inside(neighbor))
                continue;
            Position3 pos(space.coordinate2position(space.private2coord(neighbor)));
            Position3 vec((pos-center)/voxel_radius/2);
            Real r_ratio(length(pos-center)/voxel_radius/2);
            BOOST_ASSERT(r_ratio < 1.0001);
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
            const Integer nrnd((col&1)==1?8:9);
            std::pair<LatticeSpace::private_coordinate_type, bool> retval(
                    space.move_to_neighbor(private_coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.private_coord2global(retval.first).layer,
                    layer_size-1);
        }
    fout = H5::H5File("periodic_layer_1.h5", H5F_ACC_TRUNC);
    group = fout.createGroup("LatticeSpace");
    space.save(&group);
    // from layer_size-1 to 0
    layer = layer_size - 1;
    for (int row(0); row < row_size; ++row)
        for (int col(0); col < col_size; ++col)
        {
            const LatticeSpace::private_coordinate_type private_coord(
                    space.global2private_coord(Global(col, row, layer)));
            const Integer nrnd((col&1)==1?10:11);
            std::pair<LatticeSpace::private_coordinate_type, bool> retval(
                    space.move_to_neighbor(private_coord, nrnd));
            BOOST_CHECK(retval.second);
            BOOST_CHECK_EQUAL(space.private_coord2global(retval.first).layer, 0);
        }
    fout = H5::H5File("periodic_layer_2.h5", H5F_ACC_TRUNC);
    group = fout.createGroup("LatticeSpace");
    space.save(&group);
}

BOOST_AUTO_TEST_CASE(LatticeSpace_test_coordinates2)
{
    const Global g1(4, 4, 4);
    const LatticeSpace::coordinate_type c1(space.global2coord(g1));
    const LatticeSpace::private_coordinate_type pc1(space.global2private_coord(g1));
    const Global g2(space.coord2global(c1));
    const Global g3(space.private_coord2global(pc1));

    BOOST_CHECK_EQUAL(space.private2coord(space.coord2private(c1)), c1);
    BOOST_CHECK_EQUAL(space.coord2private(c1), pc1);
    BOOST_CHECK_EQUAL(space.private2coord(pc1), c1);

    BOOST_CHECK(g1.col == g2.col && g1.row == g2.row && g1.layer == g2.layer);
    BOOST_CHECK(g1.col == g3.col && g1.row == g3.row && g1.layer == g3.layer);

    const Position3 p1(space.global2position(g1));
    const Global g4(space.position2global(p1));

    BOOST_CHECK(g1.col == g4.col && g1.row == g4.row && g1.layer == g4.layer);
    BOOST_CHECK_EQUAL(c1, space.position2coordinate(p1));
    BOOST_CHECK_EQUAL(pc1, space.position2private(p1));

    const Position3 p2(space.coordinate2position(c1));
    BOOST_CHECK_EQUAL(c1, space.position2coordinate(p2));
}

BOOST_AUTO_TEST_SUITE_END()

