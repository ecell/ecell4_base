#include <ecell4/core/Model.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/STLFileIO.hpp>
#include <ecell4/sgfrd/SGFRDWorld.hpp>
#include <ecell4/sgfrd/SGFRDSimulator.hpp>

#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iomanip>

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
           << world->num_molecules(sp3) << std::endl;
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
        outstr << "t = " << std::setprecision(17) << iter->second.t()
               << ", reactant = { ";
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
    typedef ecell4::Polygon polygon_type;
    typedef ecell4::sgfrd::SGFRDWorld     world_type;
    typedef ecell4::sgfrd::SGFRDSimulator simulator_type;

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

    if(key_missing(input, "system_size")) return 1;
    const ecell4::Real     L(boost::lexical_cast<ecell4::Real>(input["system_size"]));
    const ecell4::Real3    edge_lengths(L, L, L);
    const ecell4::Integer3 matrix_sizes(3, 3, 3);
    const ecell4::Real     volume(L * L * L);

    if(key_missing(input, "polygon")){return 1;}
    boost::shared_ptr<polygon_type> polygon = boost::make_shared<polygon_type>(
        polygon_type(edge_lengths, ecell4::read_stl_format(input["polygon"], ecell4::STLFormat::Ascii)));

    boost::shared_ptr<ecell4::NetworkModel>
        model(new ecell4::NetworkModel());

    if(key_missing(input, "DA")){return 1;}
    if(key_missing(input, "RA")){return 1;}
    ecell4::Species sp1(std::string("A"), boost::lexical_cast<ecell4::Real>(input["RA"]), boost::lexical_cast<ecell4::Real>(input["DA"]));
    model->add_species_attribute(sp1);
    if(key_missing(input, "DB")){return 1;}
    if(key_missing(input, "RB")){return 1;}
    ecell4::Species sp2(std::string("B"), boost::lexical_cast<ecell4::Real>(input["RB"]), boost::lexical_cast<ecell4::Real>(input["DB"]));
    model->add_species_attribute(sp2);
    if(key_missing(input, "DC")){return 1;}
    if(key_missing(input, "RC")){return 1;}
    ecell4::Species sp3(std::string("C"), boost::lexical_cast<ecell4::Real>(input["RC"]), boost::lexical_cast<ecell4::Real>(input["DC"]));
    model->add_species_attribute(sp3);

    if(key_missing(input, "k_bind"))  {return 1;}
    if(key_missing(input, "k_unbind")){return 1;}
    const ecell4::Real k_bind   = boost::lexical_cast<ecell4::Real>(input["k_bind"]),
                       k_unbind = boost::lexical_cast<ecell4::Real>(input["k_unbind"]);

    model->add_reaction_rule(ecell4::create_binding_reaction_rule(sp2, sp3, sp1, k_bind));
    model->add_reaction_rule(ecell4::create_unbinding_reaction_rule(sp1, sp2, sp3, k_unbind));

    boost::shared_ptr<ecell4::RandomNumberGenerator> rng =
        boost::make_shared<ecell4::GSLRandomNumberGenerator>();
    if(key_missing(input, "seed")){return 1;}
    rng->seed(boost::lexical_cast<unsigned int>(input["seed"]));

    boost::shared_ptr<world_type> world =
        boost::make_shared<world_type>(edge_lengths, matrix_sizes, rng, polygon);

    if(key_missing(input, "num_A")){return 1;}
    if(key_missing(input, "num_B")){return 1;}
    if(key_missing(input, "num_C")){return 1;}
    world->add_molecules(sp1, boost::lexical_cast<std::size_t>(input["num_A"]));
    world->add_molecules(sp2, boost::lexical_cast<std::size_t>(input["num_B"]));
    world->add_molecules(sp3, boost::lexical_cast<std::size_t>(input["num_C"]));

    if(key_missing(input, "trajectory")){return 1;}
    if(key_missing(input, "species"))   {return 1;}
    if(key_missing(input, "reaction"))  {return 1;}
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

    if(key_missing(input, "bd_dt"))          {return 1;}
    if(key_missing(input, "reaction_length")){return 1;}
    if(key_missing(input, "log_file"))       {return 1;}
    simulator_type sim(world, model,
        boost::lexical_cast<ecell4::Real>(input["bd_dt"]),
        boost::lexical_cast<ecell4::Real>(input["reaction_length"]),
        input["log_file"]);
    SGFRD_TRACE(std::cerr << "log_file = " << input["log_file"] << std::endl;)
    sim.initialize();

    snapshot_output(traj, world);
    species_output(spec, world, sp1, sp2, sp3);

    if(key_missing(input, "output_dt")){return 1;}
    if(key_missing(input, "gfrd_step")){return 1;}
    const ecell4::Real dt       = boost::lexical_cast<ecell4::Real>(input["output_dt"]);
    const std::size_t  num_step = boost::lexical_cast<std::size_t>(input["gfrd_step"]);
    for(std::size_t i=0; i<num_step; ++i)
    {
        while(sim.next_event_time() <= i * dt){sim.step();}
        assert(sim.next_event_time() > i * dt);
        sim.finalize(i * dt);

        snapshot_output(traj, world);
        species_output(spec, world, sp1, sp2, sp3);
        reac.close();
        reac.open(input["reaction"].c_str(), std::ios::trunc);
        reaction_output(reac, sim);

        sim.initialize();
    }
    sim.finalize();
    snapshot_output(traj, world);
    species_output(spec, world, sp1, sp2, sp3);

    reac.close();
    reac.open(input["reaction"].c_str(), std::ios::trunc);
    reaction_output(reac, sim);

    std::cerr << "STAT: reason for forming multi:\n";
    std::cerr << "Single Concical Shell Failed: "
              << sim.stat_multi_reason.show_percent(ecell4::sgfrd::SingleConicalFailed) << '\n';
    std::cerr << "Pair Shell Failed           : "
              << sim.stat_multi_reason.show_percent(ecell4::sgfrd::PairFailed) << '\n';
    std::cerr << '\n';
    std::cerr << "STAT: fired events: total = " << sim.stat_fired_events.total() << '\n';
    std::cerr << "Single Circular: "
              << sim.stat_fired_events.show_percent(ecell4::sgfrd::FireSingleCircular) << '\n';
    std::cerr << "Single Conical : "
              << sim.stat_fired_events.show_percent(ecell4::sgfrd::FireSingleConical) << '\n';
    std::cerr << "Pair           : "
              << sim.stat_fired_events.show_percent(ecell4::sgfrd::FirePair) << '\n';
    std::cerr << "Multi          : "
              << sim.stat_fired_events.show_percent(ecell4::sgfrd::FireMulti) << '\n';
    std::cerr << "Birth          : "
              << sim.stat_fired_events.show_percent(ecell4::sgfrd::FireBirth) << '\n';
    std::cerr << '\n';
    std::cerr << "STAT: formed multi size:\n";
    const std::vector<std::size_t> ks = sim.stat_multi_size.list_all_kinds();
    for(typename std::vector<std::size_t>::const_iterator
        i(ks.begin()), e(ks.end()); i!=e; ++i)
    {
        std::cerr << std::setw(8) << *i << " : "
                  << sim.stat_multi_size.show_percent(*i) << '\n';
    }
    std::cerr << std::flush;

    return 0;
}
