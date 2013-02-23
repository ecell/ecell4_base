#define BOOST_TEST_MODULE "SpatiocyteSimulator_test"
#define BOOST_TEST_NO_LIB

#include <boost/test/included/unit_test.hpp>

#include <ecell4/core/Position3.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include "../SpatiocyteSimulator.hpp"

using namespace ecell4;
using namespace ecell4::spatiocyte;


BOOST_AUTO_TEST_CASE(SpatiocyteSimulator_test_constructor)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(1e-8);
    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-9");
    const Real k1(0.0), k2(0.0);
    ecell4::Species sp1("A", D, radius), sp2("B", D, radius), sp3("C", D, radius);
    ReactionRule rr1(create_binding_reaction_rule(sp1, sp2, sp3, k1)),
        rr2(create_unbinding_reaction_rule(sp3, sp1, sp2, k2));

    boost::shared_ptr<ecell4::NetworkModel> model(new ecell4::NetworkModel());
    model->add_species(sp1);
    model->add_species(sp2);
    model->add_species(sp3);
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);

    boost::shared_ptr<SpatiocyteWorld> world(
        new SpatiocyteWorld(edge_lengths, voxel_radius));
    world->add_species(sp1);
    world->add_species(sp2);
    world->add_species(sp3);
    world->add_molecules(sp1, N / 2);
    world->add_molecules(sp2, N / 2);

    SpatiocyteSimulator target(model, world);
}

BOOST_AUTO_TEST_CASE(SpatiocyteSimulator_test_step)
{
    const Real L(1e-6);
    const Position3 edge_lengths(L, L, L);
    const Real voxel_radius(1e-8);
    const Integer N(60);

    const std::string D("1e-12"), radius("2.5e-9");
    const Real k1(0.0), k2(0.0);
    ecell4::Species sp1("A", D, radius), sp2("B", D, radius), sp3("C", D, radius);
    ReactionRule rr1(create_binding_reaction_rule(sp1, sp2, sp3, k1)),
        rr2(create_unbinding_reaction_rule(sp3, sp1, sp2, k2));

    boost::shared_ptr<ecell4::NetworkModel> model(new ecell4::NetworkModel());
    model->add_species(sp1);
    model->add_species(sp2);
    model->add_species(sp3);
    model->add_reaction_rule(rr1);
    model->add_reaction_rule(rr2);

    boost::shared_ptr<SpatiocyteWorld> world(
        new SpatiocyteWorld(edge_lengths, voxel_radius));
    world->add_species(sp1);
    world->add_species(sp2);
    world->add_species(sp3);
    world->add_molecules(sp1, N / 2);
    world->add_molecules(sp2, N / 2);

    SpatiocyteSimulator target(model, world);

    world->add_molecules(sp1, N / 2);

    BOOST_ASSERT(world->num_molecules(sp1) == N);
    BOOST_ASSERT(world->num_molecules(sp2) == N / 2);

    target.step();
}
