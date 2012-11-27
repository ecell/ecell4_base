#include <string>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include "../../BDSimulator.cpp"

using namespace ecell4;
using namespace ecell4::bd;


int main(int argc, char** argv)
{
    const Real L(1e-6);
    std::string D("5e-12"), radius("5e-9");
    Position3 const edge_lengths(L, L, L);

    boost::shared_ptr<Model> model(new NetworkModel());
    Species sp1("A");
    sp1.set_attribute("D", D);
    sp1.set_attribute("radius", radius);
    (*model).add_species(sp1);

    boost::shared_ptr<BDWorld> world(new BDWorld(edge_lengths));
    Particle const p1(
        sp1, Position3(0, 0, 0),
        std::atof(radius.c_str()), std::atof(D.c_str()));
    (*world).new_particle(p1);

    GSLRandomNumberGenerator rng;

    BDSimulator target(model, world, rng);
    target.step();
}
