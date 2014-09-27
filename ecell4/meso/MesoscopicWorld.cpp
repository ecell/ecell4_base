#include <ecell4/core/SerialIDGenerator.hpp>
#include "MesoscopicWorld.hpp"

namespace ecell4
{

namespace meso
{

MoleculeInfo MesoscopicWorld::get_molecule_info(const Species& sp) const
{
    const bool with_D(sp.has_attribute("D"));
    Real D(0.0);

    if (with_D)
    {
        D = std::atof(sp.get_attribute("D").c_str());
    }
    else
    {
        if (boost::shared_ptr<Model> bound_model = lock_model())
        {
            Species attributed(bound_model->apply_species_attributes(sp));
            if (attributed.has_attribute("D"))
            {
                D = std::atof(attributed.get_attribute("D").c_str());
            }
        }
    }

    MoleculeInfo info = {D};
    return info;
}

std::vector<std::pair<ParticleID, Particle> >
    MesoscopicWorld::list_particles() const
{
    SerialIDGenerator<ParticleID> pidgen;
    const std::vector<Species>& species_list(species());
    const Position3 lengths(subvolume_edge_lengths());

    std::vector<std::pair<ParticleID, Particle> > retval;
    for (std::vector<Species>::const_iterator i(species_list.begin());
        i != species_list.end(); ++i)
    {
        MoleculeInfo info(get_molecule_info(*i));
        for (coordinate_type j(0); j < num_subvolumes(); ++j)
        {
            const Integer num(num_molecules_exact(*i, j));
            const Global g(coord2global(j));

            for (Integer k(0); k < num; ++k)
            {
                const Position3 pos(
                    rng_->uniform(g.col * lengths[0], (g.col + 1) * lengths[0]),
                    rng_->uniform(g.row * lengths[1], (g.row + 1) * lengths[1]),
                    rng_->uniform(g.layer * lengths[2], (g.layer + 1) * lengths[2]));
                retval.push_back(
                    std::make_pair(pidgen(), Particle(*i, pos, 0.0, info.D)));
            }
        }
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    MesoscopicWorld::list_particles_exact(const Species& sp) const
{
    SerialIDGenerator<ParticleID> pidgen;
    const Position3 lengths(subvolume_edge_lengths());

    std::vector<std::pair<ParticleID, Particle> > retval;
    MoleculeInfo info(get_molecule_info(sp));
    for (coordinate_type j(0); j < num_subvolumes(); ++j)
    {
        const Integer num(num_molecules_exact(sp, j));
        const Global g(coord2global(j));

        for (Integer k(0); k < num; ++k)
        {
            const Position3 pos(
                rng_->uniform(g.col * lengths[0], (g.col + 1) * lengths[0]),
                rng_->uniform(g.row * lengths[1], (g.row + 1) * lengths[1]),
                rng_->uniform(g.layer * lengths[2], (g.layer + 1) * lengths[2]));
            retval.push_back(
                std::make_pair(pidgen(), Particle(sp, pos, 0.0, info.D)));
        }
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    MesoscopicWorld::list_particles(const Species& sp) const
{
    SerialIDGenerator<ParticleID> pidgen;
    const Position3 lengths(subvolume_edge_lengths());

    std::vector<std::pair<ParticleID, Particle> > retval;
    MoleculeInfo info(get_molecule_info(sp));
    for (coordinate_type j(0); j < num_subvolumes(); ++j)
    {
        const Integer num(num_molecules(sp, j));
        const Global g(coord2global(j));

        for (Integer k(0); k < num; ++k)
        {
            const Position3 pos(
                rng_->uniform(g.col * lengths[0], (g.col + 1) * lengths[0]),
                rng_->uniform(g.row * lengths[1], (g.row + 1) * lengths[1]),
                rng_->uniform(g.layer * lengths[2], (g.layer + 1) * lengths[2]));
            retval.push_back(
                std::make_pair(pidgen(), Particle(sp, pos, 0.0, info.D)));
        }
    }
    return retval;
}

const Position3& MesoscopicWorld::edge_lengths() const
{
    return cs_->edge_lengths();
}

const Position3 MesoscopicWorld::subvolume_edge_lengths() const
{
    return cs_->subvolume_edge_lengths();
}

const Real& MesoscopicWorld::t() const
{
    return cs_->t();
}

void MesoscopicWorld::set_t(const Real& t)
{
    cs_->set_t(t);
}

Real MesoscopicWorld::get_value(const Species& sp) const
{
    return cs_->get_value(sp);
}

Real MesoscopicWorld::get_value_exact(const Species& sp) const
{
    return cs_->get_value_exact(sp);
}

const Integer MesoscopicWorld::num_subvolumes() const
{
    return cs_->num_subvolumes();
}

const Real MesoscopicWorld::subvolume() const
{
    return cs_->subvolume();
}

const Real MesoscopicWorld::volume() const
{
    return cs_->volume();
}

MesoscopicWorld::coordinate_type MesoscopicWorld::global2coord(const Global& g) const
{
    return cs_->global2coord(g);
}

Global MesoscopicWorld::coord2global(const MesoscopicWorld::coordinate_type& c) const
{
    return cs_->coord2global(c);
}

Global MesoscopicWorld::position2global(const Position3& pos) const
{
    return cs_->position2global(pos);
}

Integer MesoscopicWorld::num_molecules(const Species& sp) const
{
    return cs_->num_molecules(sp);
}

Integer MesoscopicWorld::num_molecules_exact(const Species& sp) const
{
    return cs_->num_molecules_exact(sp);
}

Integer MesoscopicWorld::num_molecules(
    const Species& sp, const MesoscopicWorld::coordinate_type& c) const
{
    return cs_->num_molecules(sp, c);
}

Integer MesoscopicWorld::num_molecules_exact(
    const Species& sp, const MesoscopicWorld::coordinate_type& c) const
{
    return cs_->num_molecules_exact(sp, c);
}

void MesoscopicWorld::add_molecules(
    const Species& sp, const Integer& num, const MesoscopicWorld::coordinate_type& c)
{
    cs_->add_molecules(sp, num, c);
}

void MesoscopicWorld::remove_molecules(
    const Species& sp, const Integer& num, const MesoscopicWorld::coordinate_type& c)
{
    cs_->remove_molecules(sp, num, c);
}

const std::vector<Species>& MesoscopicWorld::species() const
{
    return cs_->species();
}

std::vector<Species> MesoscopicWorld::list_species() const
{
    return cs_->list_species();
}

} // meso

} // ecell4
