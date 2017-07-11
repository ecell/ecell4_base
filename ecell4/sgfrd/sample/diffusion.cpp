#include <ecell4/core/Model.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/STLFileReader.hpp>
#include <ecell4/core/STLPolygonAdapter.hpp>

#include <ecell4/sgfrd/polygon_traits.hpp>
#include <ecell4/sgfrd/SGFRDWorld.hpp>
#include <ecell4/sgfrd/SGFRDSimulator.hpp>

#include <boost/lexical_cast.hpp>

void trajectory_output(const std::vector<ecell4::ParticleID>& pids,
    const boost::shared_ptr<ecell4::sgfrd::SGFRDWorld>& world)
{
    for(std::vector<ecell4::ParticleID>::const_iterator
        iter = pids.begin(); iter != pids.end(); ++iter)
    {
        std::cout << world->t() << ' '
                  << world->get_particle(*iter).second.position()[0] << ' '
                  << world->get_particle(*iter).second.position()[1] << ' '
                  << world->get_particle(*iter).second.position()[2] << '\n';
    }
    std::cout << "\n\n";
    return;
}

int main(int argc, char **argv)
{
    typedef ecell4::sgfrd::polygon_traits       polygon_traits;
    typedef ecell4::Polygon<polygon_traits>     polygon_type;
    typedef ecell4::sgfrd::SGFRDWorld           world_type;
    typedef ecell4::sgfrd::SGFRDSimulator       simulator_type;
    typedef typename polygon_type::face_id_type face_id_type;

    if(argc != 4)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <polygon-file.stl> <system size> <num-gfrd-step>"
                  << std::endl;
        return 1;
    }

    const std::string stlname(argv[1]);
    ecell4::STLFileReader reader;
    ecell4::STLPolygonAdapter<polygon_traits> adapter;
    boost::shared_ptr<polygon_type> polygon =
        adapter.make_polygon(reader.read(stlname, ecell4::STLFileReader::Ascii));

    const ecell4::Real     L(boost::lexical_cast<ecell4::Real>(std::string(argv[2])));
    const ecell4::Real3    edge_lengths(L, L, L);
    const ecell4::Integer3 matrix_sizes(3, 3, 3);
    const ecell4::Real     volume(L * L * L);

    boost::shared_ptr<ecell4::NetworkModel>
        model(new ecell4::NetworkModel());

    ecell4::Species sp1(std::string("A"),
         /* radius = */ std::string("1.e-2"),
         /* D = */      std::string("1.e-2"));
    model->add_species_attribute(sp1);

    boost::shared_ptr<ecell4::RandomNumberGenerator> rng =
        boost::make_shared<ecell4::GSLRandomNumberGenerator>();
    rng->seed((unsigned long int)123456);

    boost::shared_ptr<world_type> world =
        boost::make_shared<world_type>(edge_lengths, matrix_sizes, polygon, rng);

    const std::size_t num_particle = polygon->num_triangles();

    const std::vector<face_id_type> fids = polygon->list_face_id();
    std::vector<ecell4::ParticleID> pids; pids.reserve(num_particle);

    for(std::size_t np = 0; np < num_particle; ++np)
    {
        ecell4::Particle p(sp1, ecell4::centroid(polygon->triangle_at(fids.at(np))),
                           1e-2, 1e-2);

        //TODO add world::new_particle(Species, pos)
        const std::pair<std::pair<ecell4::ParticleID, ecell4::Particle>, bool> newp =
            world->new_particle(p, fids.at(np));
        assert(newp.second);
        pids.push_back(newp.first.first);
    }

    simulator_type sim(world, model);
    sim.initialize();
    trajectory_output(pids, world);

    const ecell4::Real dt = 0.001;
    const std::size_t num_step = boost::lexical_cast<std::size_t>(std::string(argv[3]));
    for(std::size_t i=0; i<num_step; ++i)
    {
        while(sim.step(i * dt)){}
        trajectory_output(pids, world);
   }
    sim.finalize();
    trajectory_output(pids, world);

    return 0;
}
