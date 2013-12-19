#include "LatticeSimulator.hpp"

namespace ecell4
{

namespace lattice
{

void LatticeSimulator::initialize()
{
    //(*world_).initialize();
    dt_ = 0.1;
    // dt = (4 * R^2) / (6 * D) for each species
}

void LatticeSimulator::step()
{
    boost::shared_ptr<GSLRandomNumberGenerator> rng((*world_).rng());

    std::vector<Species> species((*world_).list_species());
    for (std::vector<Species>::iterator s_itr(species.begin());
            s_itr != species.end(); ++s_itr)
    {
        MolecularTypeBase* mt((*world_).get_molecular_type(*s_itr));
        //shuffle(*rng, mt->voxels());
        for (MolecularType::container_type::iterator itr(mt->begin());
                itr != mt->end(); ++itr)
        {
            Coord from_coord((*itr).first);
            Coord to_coord((*world_).get_neighbor(from_coord,
                        (*rng).uniform_int(0, 11)));
            (*world_).move(from_coord, to_coord);
            /*
            if ((*world_).get_molecular_type(to_coord)->is_vacant())
            {
                (*world_).update_particle(*itr, mt, to_coord);
            }
            */
        }
    }
    (*world_).set_t(t() + dt());
}

bool LatticeSimulator::step(const Real& upto)
{
    Integer count((Integer)(upto / dt()) - num_steps());
    for (Integer i(0); i < count; ++i)
    {
        step();
    }
}

} // lattice

} // ecell4
