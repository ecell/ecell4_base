#include <ecell4/core/NetworkModel.hpp>
#include "../LatticeSimulator.hpp"
#include <ecell4/core/Rod.hpp>

using namespace ecell4;

int main()
{
    const std::string radius("1e-8"),
                      cytoplasm("Cytoplasm"),
                      membrane("Membrane");
    Species minDadp("MinDadp", radius, "16e-12", cytoplasm),
            minDEE("MinDEE", radius, "0.02e-12", membrane);

    boost::shared_ptr<NetworkModel> model(new NetworkModel);
    model->add_species_attribute(Species(cytoplasm, radius, "0"));
    model->add_species_attribute(Species(membrane, radius, "0"));
    model->add_species_attribute(minDadp);
    model->add_species_attribute(minDEE);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<lattice::LatticeWorld> world(
            new lattice::LatticeWorld(Real3(4.6e-6, 1.1e-6, 1.1e-6), 1e-8, rng));

    boost::shared_ptr<const Rod> rod(new Rod(3.5e-6, 0.5e-6, Real3(2.3e-6, 0.55e-6, 0.55e-6)));
    boost::shared_ptr<const RodSurface> surface_ptr(
            new RodSurface(3.5e-6, 0.5e-6, Real3(2.3e-6, 0.55e-6, 0.55e-6)));

    world->add_structure(Species(cytoplasm), rod);
    world->add_structure(Species(membrane), surface_ptr);

    world->add_molecules(minDadp, 1300);
    world->add_molecules(minDEE, 700);

    lattice::LatticeSimulator sim(model, world);

    sim.run(1e-1);
}
