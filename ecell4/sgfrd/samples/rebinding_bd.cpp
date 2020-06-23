// # input file
// polygon          = plane10x10.stl # triangles as STL format
// system_size      = 10.0           # boundary size
// bd_dt            = 2.5e-7         # overall BD dt (< (0.1 * R_min)^2 / D_max)
// reaction_length  = 0.0005         # reaction length (0.1 * R_min)
//
// DA = 1.0
// RA = 0.005
// DB = 0.0
// RB = 0.005
//
// k_bind   = 1.256
// k_unbind = 0.0
//
// gfrd_step = 100
// output_dt = 0.002
// seed      = 123456789
//
// trajectory = rebinding_bd_trajectory.dat
// species    = rebinding_bd_species.dat
// reaction   = rebinding_bd_reaction_record.dat
// log_file   = rebinding_bd_trace.log
//
#include <ecell4/core/Model.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/STLFileIO.hpp>
#include <ecell4/sgfrd/SGFRDWorld.hpp>
#include <ecell4/sgfrd/SGFRDSimulator.hpp>
#include <ecell4/sgfrd/BDSimulator.hpp>

#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iomanip>

void snapshot_output(std::ofstream& outstr,
        const std::shared_ptr<ecell4::sgfrd::SGFRDWorld>& world)
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
        const std::shared_ptr<ecell4::sgfrd::SGFRDWorld>& world,
        const ecell4::Species& sp1, const ecell4::Species& sp2)
{
    outstr << world->t() << ' '
           << world->num_molecules(sp1) << ' '
           << world->num_molecules(sp2) << std::endl;
    return;
}

template<typename SimulatorT>
void reaction_output(std::ofstream& outstr, const SimulatorT& sim)
{
    typedef typename SimulatorT::reaction_rule_type reaction_rule_type;
    typedef typename SimulatorT::reaction_info_type reaction_info_type;
    typedef std::vector<std::pair<reaction_rule_type, reaction_info_type> >
            reaction_record_type;

    reaction_record_type const& rr = sim.last_reactions();
    for(typename reaction_record_type::const_iterator
            iter = rr.begin(), end_ = rr.end(); iter != end_; ++iter)
    {
        outstr << "rule: " << iter->first.as_string() << " : ";
        outstr << "t = " << std::setprecision(17) << iter->second.t()
               << ", reactant = { ";
        for(typename reaction_info_type::container_type::const_iterator
                r_iter = iter->second.reactants().begin(),
                r_end = iter->second.reactants().end(); r_iter != r_end; ++r_iter)
        {
            outstr << r_iter->first;
        }
        outstr << "}, products = { ";
        for(typename reaction_info_type::container_type::const_iterator
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
    typedef ecell4::sgfrd::BDSimulator    simulator_type;
//  typedef ecell4::sgfrd::SGFRDSimulator simulator_type;

    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " input.inp" << std::endl;
        return 1;
    }

    // -----------------------------------------------------------------
    // read input file
    const std::string inpname(argv[1]);
    std::cout << "input file = " << inpname << std::endl;

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
        if(line.empty() || line[0] == '#') {continue;}

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

    // -----------------------------------------------------------------
    // setup world, model, and simulator

    if(key_missing(input, "system_size")){return 1;}
    const ecell4::Real     L(boost::lexical_cast<ecell4::Real>(input["system_size"]));
    const ecell4::Real3    edge_lengths(L, L, L);
    const ecell4::Integer3 matrix_sizes(3, 3, 3);
    const ecell4::Real     volume(L * L * L);

    if(key_missing(input, "polygon")){return 1;}
    std::shared_ptr<polygon_type> polygon = std::make_shared<polygon_type>(
            polygon_type(edge_lengths, ecell4::read_stl_format(
                    input["polygon"], ecell4::STLFormat::Ascii))
        );

    const std::vector<ecell4::FaceID> fids = polygon->list_face_ids();

    std::shared_ptr<ecell4::NetworkModel> model(new ecell4::NetworkModel());

    // A + B -> B, irreversible

    if(key_missing(input, "DA")){return 1;}
    if(key_missing(input, "RA")){return 1;}
    const ecell4::Real RA = boost::lexical_cast<ecell4::Real>(input["RA"]);
    const ecell4::Real DA = boost::lexical_cast<ecell4::Real>(input["DA"]);
    ecell4::Species sp1(std::string("A"), RA, DA);
    model->add_species_attribute(sp1);

    if(key_missing(input, "DB")){return 1;}
    if(key_missing(input, "RB")){return 1;}
    const ecell4::Real RB = boost::lexical_cast<ecell4::Real>(input["RB"]);
    const ecell4::Real DB = boost::lexical_cast<ecell4::Real>(input["DB"]);
    ecell4::Species sp2(std::string("B"), RB, DB);
    model->add_species_attribute(sp2);

    if(key_missing(input, "k_bind"))  {return 1;}
    const ecell4::Real k_bind = boost::lexical_cast<ecell4::Real>(input["k_bind"]);
    model->add_reaction_rule(ecell4::create_binding_reaction_rule(sp1, sp2, sp2, k_bind));

    std::shared_ptr<ecell4::RandomNumberGenerator> rng =
        std::make_shared<ecell4::GSLRandomNumberGenerator>();
    if(key_missing(input, "seed")){return 1;}
    rng->seed(boost::lexical_cast<unsigned int>(input["seed"]));

    std::shared_ptr<world_type> world =
        std::make_shared<world_type>(edge_lengths, matrix_sizes, polygon, rng);

    // -----------------------------------------------------------------
    // put particle that are in contact with each other.
    const ecell4::Triangle face = polygon->triangle_at(fids.front());
    const ecell4::Real3    com  =
        (face.vertex_at(0) + face.vertex_at(1) + face.vertex_at(2)) / 3.0;

    // XXX assuming planer surface along 0,0 to 10,10
    const ecell4::Real safety = 1.0 + 1e-6;
    const ecell4::Real distA_safety = RA * safety * std::sqrt(0.5);
    const ecell4::Real distB_safety = RB * safety * std::sqrt(0.5);

    ecell4::::FaceID fid_A = fids.front();
    ecell4::::FaceID fid_B = fids.front();
    ecell4::Real3 pos_A = com + ecell4::Real3(distA_safety, distA_safety, 0);
    ecell4::Real3 pos_B = com - ecell4::Real3(distB_safety, distB_safety, 0);

    std::cout << std::setprecision(17) << ecell4::length(pos_A - pos_B) << " > " << RA+RB << std::endl;


    // move particle once to remove the effect of reaction volume.
    // here, assuming the polygon is planer surface.
    bool move_A;
    if(DA <= 0.0 && DB <= 0.0)
    {
        std::cerr << "[error] neither particle moves" << std::endl;
        return 1;
    }
    else if(DA == 0.0 && DB > 0.0)
    {
        move_A = false;
    }
    else if(DA == 0.0 && DB > 0.0)
    {
        move_A = true;
    }
    else
    {
        move_A = (rng->uniform_int(0, 1) == 0);
    }


    // read dt
    if(key_missing(input, "output_dt")){return 1;}
    const ecell4::Real dt = boost::lexical_cast<ecell4::Real>(input["output_dt"]);

    // this code here does not consider the curvature of the surface of the polygon
    // (assuming the polygon is an ideal planner surface)
    if(move_A)
    {
        const ecell4::Real  sqrt2Dt = std::sqrt(2 * DA * dt);
        const ecell4::Real3 disp(rng->gaussian(sqrt2Dt), rng->gaussian(sqrt2Dt), 0.0);

        if(length(pos_A + disp - pos_B) > (RA + RB) * safety)
        {
            std::tie(pos_A, fid_A) = ecell4::polygon::travel(
                    *polygon, std::make_pair(pos_A, fid_A), disp);
        }
    }
    else // move_B
    {
        const ecell4::Real  sqrt2Dt = std::sqrt(2 * DB * dt);
        const ecell4::Real3 disp(rng->gaussian(sqrt2Dt), rng->gaussian(sqrt2Dt), 0.0);

        if(length(pos_B + disp - pos_A) > (RA + RB) * safety)
        {
            std::tie(pos_B, fid_B) = ecell4::polygon::travel(
                    *polygon, std::make_pair(pos_B, fid_B), disp);
        }
    }

    // assign particle
    world->new_particle(ecell4::Particle(sp1, pos_A, RA, DA), fid_A);
    world->new_particle(ecell4::Particle(sp2, pos_B, RB, DB), fid_B);

    // -----------------------------------------------------------------
    // open and clear output files

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

    // -----------------------------------------------------------------
    // build simulator

    if(key_missing(input, "bd_dt"))          {return 1;}
    if(key_missing(input, "reaction_length")){return 1;}
    if(key_missing(input, "log_file"))       {return 1;}
    simulator_type sim(world, model,
        boost::lexical_cast<ecell4::Real>(input["bd_dt"]),
        boost::lexical_cast<ecell4::Real>(input["reaction_length"]),
        input["log_file"]);
    SGFRD_TRACE(std::cerr << "log_file = " << input["log_file"] << std::endl;)
    sim.initialize();

    // -----------------------------------------------------------------
    // run simulation

    if(key_missing(input, "gfrd_step")){return 1;}
    const std::size_t  num_step = boost::lexical_cast<std::size_t>(input["gfrd_step"]);
    for(std::size_t i=0; i<num_step; ++i)
    {
        while(sim.next_event_time() <= i * dt)
        {
            sim.step();
            if(world->num_molecules(sp1) == 0) {break;}
        }
        if(world->num_molecules(sp1) == 0) {break;}

        assert(sim.next_event_time() > i * dt);
        sim.finalize(i * dt);

        std::cerr << "current t = " << i*dt << std::endl;

        snapshot_output(traj, world);
        species_output(spec, world, sp1, sp2);
        reac.close();
        reac.open(input["reaction"].c_str(), std::ios::trunc);
        reaction_output(reac, sim);

        sim.initialize();
    }
    sim.finalize();
    snapshot_output(traj, world);
    species_output(spec, world, sp1, sp2);

    reac.close();
    reac.open(input["reaction"].c_str(), std::ios::trunc);
    reaction_output(reac, sim);

    return 0;
}
