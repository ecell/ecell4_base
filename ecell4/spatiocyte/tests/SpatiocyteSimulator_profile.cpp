#include <ecell4/core/NetworkModel.hpp>
#include "../SpatiocyteSimulator.hpp"
#include <ecell4/core/Rod.hpp>

using namespace ecell4;

int main()
{
    const std::string radius("1e-8"),
                      cytoplasm_str("Cytoplasm"),
                      membrane_str("Membrane");
    Species cytoplasm(cytoplasm_str, radius, "0"),
            membrane(membrane_str, radius, "0", cytoplasm_str),
            minDadp("MinDadp", radius, "16e-12", cytoplasm_str),
            minDEE("MinDEE", radius, "0.02e-12", membrane_str);

    boost::shared_ptr<NetworkModel> model(new NetworkModel);
    model->add_species_attribute(cytoplasm);
    model->add_species_attribute(membrane);
    model->add_species_attribute(minDadp);
    model->add_species_attribute(minDEE);

    boost::shared_ptr<GSLRandomNumberGenerator>
        rng(new GSLRandomNumberGenerator());
    boost::shared_ptr<spatiocyte::SpatiocyteWorld> world(
            new spatiocyte::SpatiocyteWorld(Real3(4.6e-6, 1.1e-6, 1.1e-6), 1e-8, rng));

    boost::shared_ptr<const Rod> rod(
            new Rod(3.5e-6, 0.5e-6, Real3(2.3e-6, 0.55e-6, 0.55e-6)));
    boost::shared_ptr<const RodSurface> surface(
            new RodSurface(3.5e-6, 0.5e-6, Real3(2.3e-6, 0.55e-6, 0.55e-6)));

    world->add_structure(cytoplasm, rod);
    world->add_structure(membrane, surface);

    world->add_molecules(minDadp, 1300);
    world->add_molecules(minDEE, 700);

    spatiocyte::SpatiocyteSimulator sim(model, world);

    sim.run(1e-1);
}
