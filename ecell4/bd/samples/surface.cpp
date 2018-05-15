#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/STLFileIO.hpp>

#include <ecell4/bd/BDWorld.hpp>
#include <ecell4/bd/BDSimulator.hpp>

using namespace ecell4;
using namespace ecell4::bd;

/**
 * a simple function to dump particle position(s)
 */
void print_particle_position(const BDWorld& world)
{
    const std::vector<std::pair<ParticleID, Particle> > ps(world.list_particles());
    for(std::vector<std::pair<ParticleID, Particle> >::const_iterator
            iter = ps.begin(); iter != ps.end(); ++iter)
    {
        const Real3 pos(iter->second.position());
        std::cout << iter->second.sid() << " " << std::setprecision(12)
                  << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
    }
    std::cout << std::endl << std::endl;
    return ;
}

/**
 * main function
 */
int main(int argc, char** argv)
{
    if(argc != 3)
    {
        std::cerr << "usage: " << argv[0] << " -stl=[bin|asc] [filename.stl]" << std::endl;
        return 1;
    }

    /// simulation parameters
    const Real L(10);
    std::string D("0.5"), radius("1e-2");
    const Real3 edge_lengths(L, L, L);
    const Integer3 matrix_sizes(3, 3, 3);

    std::string stltype(argv[1]);
    if(stltype.size() != 8)
    {
        std::cerr << "usage: " << argv[0] << " -stl=[bin|asc] [filename.stl]" << std::endl;
        return 1;
    }
    STLFormat::Kind kind;
    if(stltype.substr(5, 3) == "asc")
    {
        kind = STLFormat::Ascii;
    }
    else if(stltype.substr(5, 3) == "bin")
    {
        kind = STLFormat::Binary;
    }
    else
    {
        std::cerr << "usage: " << argv[0] << " -stl=[bin|asc] [filename.stl]" << std::endl;
        return 1;
    }

    /// instantiate NetworkModel
    boost::shared_ptr<NetworkModel> model(new NetworkModel());

    /// create a Species, and set its attributes
    Species sp1("A");
    sp1.set_attribute("D", D);
    sp1.set_attribute("radius", radius);

    Species sp2("B");
    sp2.set_attribute("D", D);
    sp2.set_attribute("radius", radius);

    // A -> A + A
//     ReactionRule rr1;
//     rr1.set_k(1.0);
//     rr1.add_reactant(sp1);
//     rr1.add_product(sp1);
//     rr1.add_product(sp1);
//
//     // A + A -> B
//     ReactionRule rr2;
//     rr2.set_k(1.0);
//     rr2.add_reactant(sp1);
//     rr2.add_reactant(sp1);
//     rr2.add_product(sp2);

    (*model).add_species_attribute(sp1);
//     (*model).add_species_attribute(sp2);
//     (*model).add_reaction_rule(rr1);
//     (*model).add_reaction_rule(rr2);

    boost::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());

    /// instantiate BDWorld
    const std::string stlname(argv[2]);
    boost::shared_ptr<BDWorld> world(
        new BDWorld(Polygon(edge_lengths, read_stl_format(stlname, kind)),
                    matrix_sizes, rng));
    world->bind_to(model);

    std::cout << "world instanciated" << std::endl;
    std::cout << "polygon has " << world->polygon().face_size()   << " faces\n";
    std::cout << "        and " << world->polygon().edge_size()   << " edges\n";
    std::cout << "        and " << world->polygon().vertex_size() << " vertices\n";

    /// create a Particle, and inject it into BDWorld
    std::cerr << "begin create/inject particle to BDWorld";
    for(std::size_t i=0; i<100; ++i)
    {
        BDWorld::molecule_info_type info1((*world).get_molecule_info(Species("A")));

        Polygon::FaceID fid;
        const Real3     pos = world->polygon().draw_position(rng, fid);
        Particle p(sp1, pos, info1.radius, info1.D);
        world->new_particle(p, fid);
    }
    std::cerr << "... end!" << std::endl;

    std::cerr << "begin simulator setup";
    /// instatiate BDSimulator
    BDSimulator sim(model, world);
    sim.set_dt(1e-6);
    std::cerr << "... end!" << std::endl;

    /// run and log
    std::cerr << "begin simulation..." << std::endl;
    for(unsigned int i(0); i <= 100; ++i)
    {
        while(sim.step(1e-4*i))
        {
            // do nothing!
        }
        print_particle_position(*world);
    }
    std::cerr << "end!" << std::endl;
    return 0;
}
