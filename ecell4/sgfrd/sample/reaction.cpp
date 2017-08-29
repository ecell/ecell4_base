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
#include <fstream>

void snapshot_output(std::ofstream& outstr,
        const boost::shared_ptr<ecell4::sgfrd::SGFRDWorld>& world)
{
    typedef ecell4::sgfrd::SGFRDWorld::particle_container_type container;
    container const ps = world->list_particles();
    for(container::const_iterator iter = ps.begin(); iter != ps.end(); ++iter)
    {
        outstr << world->t() << ' '
               << iter->second.species().serial() << ' '
               << iter->second.position()[0]      << ' '
               << iter->second.position()[1]      << ' '
               << iter->second.position()[2]      << '\n';
    }
    outstr << "\n\n" << std::flush;
    return;
}

void species_output(std::ofstream& outstr,
        const boost::shared_ptr<ecell4::sgfrd::SGFRDWorld>& world,
        const ecell4::Species& sp1, const ecell4::Species& sp2,
        const ecell4::Species& sp3)
{
    outstr << world->t() << ' '
           << world->num_molecules(sp1) << ' '
           << world->num_molecules(sp2) << ' '
           << world->num_molecules(sp3) << '\n';
    return;
}

void reaction_output(
        std::ofstream& outstr, const ecell4::sgfrd::SGFRDSimulator& sim)
{
    typedef ecell4::sgfrd::SGFRDSimulator::reaction_rule_type reaction_rule_type;
    typedef ecell4::sgfrd::SGFRDSimulator::reaction_info_type reaction_info_type;
    typedef std::vector<std::pair<reaction_rule_type, reaction_info_type> >
            reaction_record_type;

    reaction_record_type const& rr = sim.last_reactions();
    for(reaction_record_type::const_iterator iter = rr.begin(), end_ = rr.end();
            iter != end_; ++iter)
    {
        outstr << "rule: " << iter->first.as_string() << " : ";
        outstr << "t = " << iter->second.t() << ", reactant = { ";
        for(reaction_info_type::container_type::const_iterator
                r_iter = iter->second.reactants().begin(),
                r_end = iter->second.reactants().end(); r_iter != r_end; ++r_iter)
        {
            outstr << r_iter->first;
        }
        outstr << "}, products = { ";
        for(reaction_info_type::container_type::const_iterator
                p_iter = iter->second.products().begin(),
                p_end  = iter->second.products().end(); p_iter != p_end; ++p_iter)
        {
            outstr << p_iter->first;
        }
        outstr << '}' << std::endl;
    }

    return;
}


bool key_missing(
        std::map<std::string, std::string> const& inp, std::string const& key)
{
    if(inp.count(key) != 1)
    {
        std::cerr << "missing key " << key << " in input file" << std::endl;
        return true;
    }
    return false;
}

int main(int argc, char **argv)
{
    typedef ecell4::sgfrd::polygon_traits       polygon_traits;
    typedef ecell4::Polygon<polygon_traits>     polygon_type;
    typedef ecell4::sgfrd::SGFRDWorld           world_type;
    typedef ecell4::sgfrd::SGFRDSimulator       simulator_type;
    typedef typename polygon_type::face_id_type face_id_type;

    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " input.inp" << std::endl;
        return 1;
    }

    const std::string inpname(argv[1]);
    std::ifstream ifs(inpname.c_str());
    if(!ifs.good())
    {
        std::cerr << "file open error: " << inpname << std::endl;
        return 1;
    }

    std::map<std::string, std::string> input;
    while(!ifs.eof())
    {
        std::string line;
        std::getline(ifs, line);
        if(line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key, equal, value;
        iss >> key >> equal >> value;
        if(equal != "=")
        {
            std::cerr << "invalid line in input file: " << line << std::endl;
            return 1;
        }
        input.insert(std::make_pair(key, value));
        ifs.peek();
    }

    ecell4::STLFileReader reader;
    ecell4::STLPolygonAdapter<polygon_traits> adapter;

    if(key_missing(input, "polygon")) return 1;
    boost::shared_ptr<polygon_type> polygon =
        adapter.make_polygon(reader.read(input["polygon"], ecell4::STLFileReader::Ascii));

    if(key_missing(input, "system_size")) return 1;
    const ecell4::Real     L(boost::lexical_cast<ecell4::Real>(input["system_size"]));
    const ecell4::Real3    edge_lengths(L, L, L);
    const ecell4::Integer3 matrix_sizes(3, 3, 3);
    const ecell4::Real     volume(L * L * L);

    boost::shared_ptr<ecell4::NetworkModel>
        model(new ecell4::NetworkModel());

    if(key_missing(input, "DA")) return 1;
    if(key_missing(input, "RA")) return 1;
    ecell4::Species sp1(std::string("A"), input["RA"], input["DA"]);
    model->add_species_attribute(sp1);
    if(key_missing(input, "DB")) return 1;
    if(key_missing(input, "RB")) return 1;
    ecell4::Species sp2(std::string("B"), input["RB"], input["DB"]);
    model->add_species_attribute(sp2);
    if(key_missing(input, "DC")) return 1;
    if(key_missing(input, "RC")) return 1;
    ecell4::Species sp3(std::string("C"), input["RC"], input["DC"]);
    model->add_species_attribute(sp3);

    if(key_missing(input, "k_bind")) return 1;
    if(key_missing(input, "k_unbind")) return 1;
    const ecell4::Real k_bind   = boost::lexical_cast<ecell4::Real>(input["k_bind"]),
                       k_unbind = boost::lexical_cast<ecell4::Real>(input["k_unbind"]);

    model->add_reaction_rule(ecell4::create_binding_reaction_rule(sp2, sp3, sp1, k_bind));
    model->add_reaction_rule(ecell4::create_unbinding_reaction_rule(sp1, sp2, sp3, k_unbind));

    boost::shared_ptr<ecell4::RandomNumberGenerator> rng =
        boost::make_shared<ecell4::GSLRandomNumberGenerator>();
    rng->seed((unsigned long int)123456);

    boost::shared_ptr<world_type> world =
        boost::make_shared<world_type>(edge_lengths, matrix_sizes, polygon, rng);

    const std::size_t num_particle = polygon->num_triangles();
    world->add_molecule(sp1, num_particle);
    world->add_molecule(sp2, num_particle);
    world->add_molecule(sp3, num_particle);
    std::cout << num_particle << "*3 molecules have been added" << std::endl;

    if(key_missing(input, "trajectory")) return 1;
    if(key_missing(input, "species"))    return 1;
    if(key_missing(input, "reaction"))   return 1;
    std::ofstream traj(input["trajectory"].c_str());
    std::ofstream spec(input["species"].c_str());
    std::ofstream reac(input["reaction"].c_str());
    if(!traj.good())
    {
        std::cerr << "file open error: " << input["trajecotry"] << std::endl;
        return 1;
    }
    if(!spec.good())
    {
        std::cerr << "file open error: " << input["species"] << std::endl;
        return 1;
    }
    if(!reac.good())
    {
        std::cerr << "file open error: " << input["reaction"] << std::endl;
        return 1;
    }

    if(key_missing(input, "bd_dt")) return 1;
    if(key_missing(input, "reaction_length")) return 1;

    simulator_type sim(world, model,
        boost::lexical_cast<ecell4::Real>(input["bd_dt"]),
        boost::lexical_cast<ecell4::Real>(input["reaction_length"]));
    sim.initialize();

    snapshot_output(traj, world);
    species_output(spec, world, sp1, sp2, sp3);

    if(key_missing(input, "output_dt")) return 1;
    if(key_missing(input, "gfrd_step")) return 1;
    const ecell4::Real dt       = boost::lexical_cast<ecell4::Real>(input["output_dt"]);
    const std::size_t  num_step = boost::lexical_cast<std::size_t>(input["gfrd_step"]);
    for(std::size_t i=0; i<num_step; ++i)
    {
        while(sim.step(i * dt)){}
        snapshot_output(traj, world);
        species_output(spec, world, sp1, sp2, sp3);
    }
    sim.finalize();
    snapshot_output(traj, world);
    species_output(spec, world, sp1, sp2, sp3);
    reaction_output(reac, sim);

    return 0;
}
