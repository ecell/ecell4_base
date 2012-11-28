#include <string>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include "../../BDSimulator.cpp"

using namespace ecell4;
using namespace ecell4::bd;


void print_particle_position(BDWorld const& world, ParticleID const& pid)
{
    Position3 const pos(world.get_particle(pid).second.position());
    std::cout << std::setprecision(12) << world.t() << " : " << pos << std::endl;
}

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
    SpeciesInfo info((*world).get_species_info(sp1));
    Particle const p1(
        sp1, Position3(0, 0, 0), info.radius, info.D);
    ParticleID const pid1((*world).new_particle(p1));

    GSLRandomNumberGenerator rng;

    BDSimulator sim(model, world, rng);
    sim.set_dt(1e-6);

    const Real upto(1e-3);
    for (unsigned int i(0); i <= 10; ++i)
    {
        while (sim.step(1e-3 * i))
        {
            ; // do nothing
        }

        print_particle_position(*world, pid1);
    }
}
