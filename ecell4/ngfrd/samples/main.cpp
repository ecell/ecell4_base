#include <ecell4/ngfrd/NGFRDWorld.hpp>
#include <ecell4/ngfrd/NGFRDSimulator.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Species.hpp>

void snapshot_output(std::ofstream& ofs,
        const std::shared_ptr<ecell4::ngfrd::NGFRDWorld>& world)
{
    ofs << world->num_particles() << '\n';
    ofs << world->t() << '\n';
    for(const auto& p : world->list_particles())
    {
        ofs << p.second.species().serial() << ' '
            << p.second.position()[0]      << ' '
            << p.second.position()[1]      << ' '
            << p.second.position()[2]      << '\n';
    }
    ofs << std::flush;
    return;
}

int main()
{
    const ecell4::Real     L(10.0);
    const ecell4::Real3    edge_lengths(L, L, L);
    const ecell4::Integer3 matrix_sizes(10, 10, 10);
    const ecell4::Real     volume(L * L * L);

    std::shared_ptr<ecell4::NetworkModel> model(new ecell4::NetworkModel());
    ecell4::Species sp1(std::string("A"), /*r = */ 1e-2, /*D = */ 1e-2);
    ecell4::Species sp2(std::string("B"), /*r = */ 1e-2, /*D = */ 1e-5);
    model->add_species_attribute(sp1);
    model->add_species_attribute(sp2);

    std::shared_ptr<ecell4::RandomNumberGenerator> rng =
        std::make_shared<ecell4::GSLRandomNumberGenerator>();
    rng->seed(123456ul);

    std::shared_ptr<ecell4::ngfrd::NGFRDWorld> world = std::make_shared<ecell4::ngfrd::NGFRDWorld>(edge_lengths, matrix_sizes);
    world->add_molecules_3D(sp1, 10);
    world->add_molecules_2D(sp2, 10);

    std::ofstream traj("traj.xyz");
    snapshot_output(traj, world);

    ecell4::ngfrd::NGFRDSimulator sim(world, model);

    const ecell4::Real dt(1e-2);
    for(std::size_t i=1; i<100; ++i)
    {
        while(sim.next_event_time() <= i * dt)
        {
            std::cout << "t = " << world->t() << std::endl;
            sim.step();
        }
        assert(sim.next_event_time() > i * dt);
        sim.finalize();
        snapshot_output(traj, world);
        sim.initialize();
    }

    return 0;
}
