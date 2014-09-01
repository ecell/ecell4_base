#include "MesoscopicWorld.hpp"

namespace ecell4
{

namespace meso
{

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

std::vector<Species> MesoscopicWorld::list_species() const
{
    return cs_->list_species();
}

} // meso

} // ecell4
