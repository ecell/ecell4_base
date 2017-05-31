#include <ecell4/core/SerialIDGenerator.hpp>
#include "MesoscopicWorld.hpp"

#ifdef WIN32_MSC
#include <boost/numeric/interval/detail/msvc_rounding_control.hpp>
#endif

namespace ecell4
{

namespace meso
{

#ifdef WIN32_MSC
double round(const double x)
{
    return floor(x + 0.5);
}
#endif

MesoscopicWorld::MesoscopicWorld(const Real3& edge_lengths, const Real subvolume_length)
    : cs_(new SubvolumeSpaceVectorImpl(edge_lengths, Integer3(round(edge_lengths[0] / subvolume_length), round(edge_lengths[1] / subvolume_length), round(edge_lengths[2] / subvolume_length))))
{
    rng_ = boost::shared_ptr<RandomNumberGenerator>(
        new GSLRandomNumberGenerator());
    (*rng_).seed();
}

MesoscopicWorld::MesoscopicWorld(
    const Real3& edge_lengths, const Real subvolume_length,
    boost::shared_ptr<RandomNumberGenerator> rng)
    : cs_(new SubvolumeSpaceVectorImpl(edge_lengths, Integer3(round(edge_lengths[0] / subvolume_length), round(edge_lengths[1] / subvolume_length), round(edge_lengths[2] / subvolume_length)))), rng_(rng)
{
    ;
}

MoleculeInfo MesoscopicWorld::get_molecule_info(const Species& sp) const
{
    const bool with_D(sp.has_attribute("D"));
    const bool with_loc(sp.has_attribute("location"));

    Real D(0.0);
    std::string loc("");

    if (with_loc)
    {
        loc = sp.get_attribute_as<std::string>("location");
    }

    if (with_D)
    {
        D = sp.get_attribute_as<Real>("D");
    }
    else
    {
        if (boost::shared_ptr<Model> bound_model = lock_model())
        {
            Species newsp(bound_model->apply_species_attributes(sp));
            if (newsp.has_attribute("D"))
            {
                D = newsp.get_attribute_as<Real>("D");
            }

            if (!with_loc && newsp.has_attribute("location"))
            {
                loc = newsp.get_attribute_as<std::string>("location");
            }
        }
    }

    MoleculeInfo info = {D, loc};
    return info;
}

std::vector<std::pair<ParticleID, Particle> >
    MesoscopicWorld::list_particles() const
{
    SerialIDGenerator<ParticleID> pidgen;
    const std::vector<Species>& species_list(species());
    const Real3 lengths(subvolume_edge_lengths());

    std::vector<std::pair<ParticleID, Particle> > retval;
    for (std::vector<Species>::const_iterator i(species_list.begin());
        i != species_list.end(); ++i)
    {
        const boost::shared_ptr<PoolBase>& pool = get_pool(*i);
        for (coordinate_type j(0); j < num_subvolumes(); ++j)
        {
            const Integer num(pool->num_molecules(j));
            const Integer3 g(coord2global(j));

            for (Integer k(0); k < num; ++k)
            {
                const Real3 pos(
                    rng_->uniform(g.col * lengths[0], (g.col + 1) * lengths[0]),
                    rng_->uniform(g.row * lengths[1], (g.row + 1) * lengths[1]),
                    rng_->uniform(g.layer * lengths[2], (g.layer + 1) * lengths[2]));
                retval.push_back(
                    std::make_pair(pidgen(), Particle(*i, pos, 0.0, pool->D())));
            }
        }
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    MesoscopicWorld::list_particles_exact(const Species& sp) const
{
    SerialIDGenerator<ParticleID> pidgen;
    const Real3 lengths(subvolume_edge_lengths());

    std::vector<std::pair<ParticleID, Particle> > retval;
    if (has_species(sp))
    {
        const boost::shared_ptr<PoolBase>& pool = get_pool(sp);
        for (coordinate_type j(0); j < num_subvolumes(); ++j)
        {
            const Integer num(pool->num_molecules(j));
            const Integer3 g(coord2global(j));

            for (Integer k(0); k < num; ++k)
            {
                const Real3 pos(
                    rng_->uniform(g.col * lengths[0], (g.col + 1) * lengths[0]),
                    rng_->uniform(g.row * lengths[1], (g.row + 1) * lengths[1]),
                    rng_->uniform(g.layer * lengths[2], (g.layer + 1) * lengths[2]));
                retval.push_back(
                    std::make_pair(pidgen(), Particle(sp, pos, 0.0, pool->D())));
            }
        }
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    MesoscopicWorld::list_particles(const Species& sp) const
{
    SerialIDGenerator<ParticleID> pidgen;
    const std::vector<Species>& species_list(species());
    const Real3 lengths(subvolume_edge_lengths());

    std::vector<std::pair<ParticleID, Particle> > retval;
    // MoleculeInfo info(get_molecule_info(sp));
    // for (coordinate_type j(0); j < num_subvolumes(); ++j)
    // {
    //     const Integer num(num_molecules(sp, j));
    //     const Integer3 g(coord2global(j));

    //     for (Integer k(0); k < num; ++k)
    //     {
    //         const Real3 pos(
    //             rng_->uniform(g.col * lengths[0], (g.col + 1) * lengths[0]),
    //             rng_->uniform(g.row * lengths[1], (g.row + 1) * lengths[1]),
    //             rng_->uniform(g.layer * lengths[2], (g.layer + 1) * lengths[2]));
    //         retval.push_back(
    //             std::make_pair(pidgen(), Particle(sp, pos, 0.0, info.D)));
    //     }
    // }
    for (std::vector<Species>::const_iterator i(species_list.begin());
        i != species_list.end(); ++i)
    {
        const Integer coef(sp.count(*i));
        if (coef == 0)
        {
            continue;
        }

        const boost::shared_ptr<PoolBase>& pool = get_pool(*i);
        for (coordinate_type j(0); j < num_subvolumes(); ++j)
        {
            const Integer num(coef * pool->num_molecules(j));
            const Integer3 g(coord2global(j));

            for (Integer k(0); k < num; ++k)
            {
                const Real3 pos(
                    rng_->uniform(g.col * lengths[0], (g.col + 1) * lengths[0]),
                    rng_->uniform(g.row * lengths[1], (g.row + 1) * lengths[1]),
                    rng_->uniform(g.layer * lengths[2], (g.layer + 1) * lengths[2]));
                retval.push_back(
                    std::make_pair(pidgen(), Particle(*i, pos, 0.0, pool->D())));
            }
        }
    }
    return retval;
}

const Real3& MesoscopicWorld::edge_lengths() const
{
    return cs_->edge_lengths();
}

const Real3 MesoscopicWorld::subvolume_edge_lengths() const
{
    return cs_->subvolume_edge_lengths();
}

const Real MesoscopicWorld::t() const
{
    return cs_->t();
}

void MesoscopicWorld::set_t(const Real& t)
{
    cs_->set_t(t);
}

void MesoscopicWorld::set_value(const Species& sp, const Real value)
{
    const Integer num1 = static_cast<Integer>(value);
    const Integer num2 = num_molecules_exact(sp);
    if (num1 > num2)
    {
        add_molecules(sp, num1 - num2);
    }
    else if (num1 < num2)
    {
        remove_molecules(sp, num2 - num1);
    }
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

MesoscopicWorld::coordinate_type MesoscopicWorld::global2coord(const Integer3& g) const
{
    return cs_->global2coord(g);
}

Integer3 MesoscopicWorld::coord2global(const MesoscopicWorld::coordinate_type& c) const
{
    return cs_->coord2global(c);
}

Integer3 MesoscopicWorld::position2global(const Real3& pos) const
{
    return cs_->position2global(pos);
}

Integer MesoscopicWorld::position2coordinate(const Real3& pos) const
{
    return cs_->position2coordinate(pos);
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
    if (!cs_->has_species(sp))
    {
        this->reserve_pool(sp);
    }
    cs_->add_molecules(sp, num, c);
}

void MesoscopicWorld::remove_molecules(
    const Species& sp, const Integer& num, const MesoscopicWorld::coordinate_type& c)
{
    cs_->remove_molecules(sp, num, c);
}

std::vector<MesoscopicWorld::coordinate_type>
MesoscopicWorld::list_coordinates(const Species& sp) const
{
    return cs_->list_coordinates(sp);  // std::move?
}

std::vector<MesoscopicWorld::coordinate_type>
MesoscopicWorld::list_coordinates_exact(const Species& sp) const
{
    return cs_->list_coordinates_exact(sp);  // std::move?
}

const std::vector<Species>& MesoscopicWorld::species() const
{
    return cs_->species();
}

std::vector<Species> MesoscopicWorld::list_species() const
{
    return cs_->list_species();
}

void MesoscopicWorld::add_structure(
    const Species& sp, const boost::shared_ptr<const Shape>& shape)
{
    cs_->add_structure(sp, shape);
}

bool MesoscopicWorld::on_structure(
    const Species& sp, const coordinate_type& coord) const
{
    if (has_species(sp))
    {
        const boost::shared_ptr<PoolBase>& pool = get_pool(sp);
        return (pool->loc() == "" || cs_->check_structure(pool->loc(), coord));
    }

    const molecule_info_type minfo(get_molecule_info(sp));
    return (minfo.loc == "" || cs_->check_structure(minfo.loc, coord));
}

Real MesoscopicWorld::get_volume(const Species& sp) const
{
    return cs_->get_volume(sp);
}

} // meso

} // ecell4
