#include <string>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include "../../BDSimulator.cpp"

using namespace ecell4;
using namespace ecell4::bd;

/**
 * a simple function to dump particle position(s)
 */
void print_particle_position(BDWorld const& world, ParticleID const& pid)
{
    Position3 const pos(world.get_particle(pid).second.position());
    std::cout << std::setprecision(12) << world.t() << " : " << pos << std::endl;
}

/**
 * main function
 */
int main(int argc, char** argv)
{
    /// simulation parameters
    const Real L(1e-6);
    std::string D("5e-12"), radius("5e-9");
    Position3 const edge_lengths(L, L, L);

    /// instantiate RandomNumberGenerator
    GSLRandomNumberGenerator rng;

    /// instantiate NetworkModel
    boost::shared_ptr<Model> model(new NetworkModel());

    /// create a Species, and set its attributes
    Species sp1("A");
    sp1.set_attribute("D", D);
    sp1.set_attribute("radius", radius);
    (*model).add_species(sp1);

    Species sp2("B");
    sp2.set_attribute("D", D);
    sp2.set_attribute("radius", radius);
    (*model).add_species(sp2);

    /// instantiate BDWorld
    boost::shared_ptr<BDWorld> world(new BDWorld(edge_lengths));

    /// create a Particle, and inject it into BDWorld
    ParticleInfo info1((*world).get_particle_info(sp1));
    ParticleInfo info2((*world).get_particle_info(sp2));

    Particle const p1(
        sp1, Position3(0, 0, 0), info1.radius, info1.D);
    ParticleID const pid1((*world).new_particle(p1));

    Particle const p2(
        sp1, Position3(L * 0.5, 0, 0), info1.radius, info1.D);
    ParticleID const pid2((*world).new_particle(p2));

    Particle const p3(
        sp1, Position3(0, L * 0.5, 0), info1.radius, info1.D);
    ParticleID const pid3((*world).new_particle(p3));

    Particle const p4(
        sp2, Position3(0, 0, L * 0.5), info2.radius, info2.D);
    ParticleID const pid4((*world).new_particle(p4));

    Particle const p5(
        sp2, Position3(L * 0.1, 0, 0), info2.radius, info2.D);
    ParticleID const pid5((*world).new_particle(p5));

    Particle const p6(
        sp2, Position3(0, L * 0.1, 0), info2.radius, info2.D);
    ParticleID const pid6((*world).new_particle(p6));

    Particle const p7(
        sp2, Position3(0, 0, L * 0.1), info2.radius, info2.D);
    ParticleID const pid7((*world).new_particle(p7));

    /// instatiate BDSimulator
    BDSimulator sim(model, world, rng);
    sim.set_dt(1e-6);

    /// run and log by the millisecond
    for (unsigned int i(0); i <= 10; ++i)
    {
        while (sim.step(1e-3 * i))
        {
            ; // do nothing
        }

        print_particle_position(*world, pid1);
    }
    std::string hoge("hoge.h5");
    sim.save_space(hoge);
}
