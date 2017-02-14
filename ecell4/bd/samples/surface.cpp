#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/NetworkModel.hpp>

#include <ecell4/bd/BDPolygon.hpp>
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
    /// simulation parameters
    const Real L(1e-6);
    std::string D("5e0"), radius("1e-1");
    const Real3 edge_lengths(L, L, L);
    const Integer3 matrix_sizes(3, 3, 3);

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
    ReactionRule rr1;
    rr1.set_k(1.);
    rr1.add_reactant(sp1);
    rr1.add_product(sp1);
    rr1.add_product(sp1);

    // A + A -> B
    ReactionRule rr2;
    rr2.set_k(1.);
    rr2.add_reactant(sp1);
    rr2.add_reactant(sp1);
    rr2.add_product(sp2);

    (*model).add_species_attribute(sp1);
    (*model).add_reaction_rule(rr1);
    (*model).add_reaction_rule(rr2);

    boost::shared_ptr<RandomNumberGenerator> rng(new GSLRandomNumberGenerator());

    /// instantiate BDWorld
    boost::shared_ptr<BDWorld> world(new BDWorld(edge_lengths, matrix_sizes, rng));
    world->bind_to(model);

    // create polygon
    std::ifstream poly("polygon.xyz");
    if(!poly.good())
    {
        std::cout << "file open error" << std::endl;
        return 1;
    }

    std::cerr << "begin reading polygon.xyz";
    BDPolygon polygon;
    std::size_t num_faces=0;
    while(!poly.eof())
    {
        std::string line;
        Real3 v1, v2, v3;
        {
            Real x, y, z;
            std::getline(poly, line);
            std::istringstream iss(line);
            iss >> x >> y >> z;
            v1 = Real3(x, y, z);
        }
        {
            Real x, y, z;
            std::getline(poly, line);
            std::istringstream iss(line);
            iss >> x >> y >> z;
            v2 = Real3(x, y, z);
        }
        {
            Real x, y, z;
            std::getline(poly, line);
            std::istringstream iss(line);
            iss >> x >> y >> z;
            v3 = Real3(x, y, z);
        }

        Triangle f(v1, v2, v3);
        polygon.add_face(f);
        poly.peek();
    }
    std::cerr << "... end!" << std::endl;

    std::cerr << "begin polygon initialization";
    world->set_polygon(polygon);
    std::cerr << "... end!" << std::endl;

    const std::size_t fid=0;
    Triangle f = world->get_face(fid);
    const Real3 initial = (f.vertex_at(0) + f.vertex_at(1) + f.vertex_at(2)) / 3.;
    std::cerr << "initial position = " << initial << std::endl;

    /// create a Particle, and inject it into BDWorld
    std::cerr << "begin create/inject particle to BDWorld";
    BDWorld::molecule_info_type info1((*world).get_molecule_info(Species("A")));
    const Particle p1(sp1, initial, info1.radius, info1.D);
    const ParticleID pid1((*world).new_particle(p1, fid).first.first);
//     world->save("test_bd.h5");
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
        while(sim.step(5e-2*i))
        {
            // do nothing
        }
        print_particle_position(*world);
    }
    std::cerr << "end!" << std::endl;
    return 0;
}
