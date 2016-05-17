#include <stdexcept>
#include <fstream>

#include "SpatiocyteWorld.hpp"

namespace ecell4
{

namespace spatiocyte
{

SpatiocyteWorld* create_spatiocyte_world_cell_list_impl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const Integer3& matrix_sizes,
    const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return new SpatiocyteWorld(
        new LatticeSpaceCellListImpl(edge_lengths, voxel_radius, matrix_sizes), rng);
}

SpatiocyteWorld* create_spatiocyte_world_vector_impl(
    const Real3& edge_lengths, const Real& voxel_radius,
    const boost::shared_ptr<RandomNumberGenerator>& rng)
{
    return new SpatiocyteWorld(
        new LatticeSpaceVectorImpl(edge_lengths, voxel_radius), rng);
}

const Real SpatiocyteWorld::t() const
{
    return (*space_).t();
}

void SpatiocyteWorld::set_t(const Real& t)
{
    (*space_).set_t(t);
}

const Real3& SpatiocyteWorld::edge_lengths() const
{
    return (*space_).edge_lengths();
}

const Real SpatiocyteWorld::volume() const
{
    return (*space_).volume();
}

Integer SpatiocyteWorld::num_species() const
{
    return (*space_).num_species();
}

void SpatiocyteWorld::set_value(const Species& sp, const Real value)
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

bool SpatiocyteWorld::has_species(const Species &sp) const
{
    return (*space_).has_species(sp);
}

// bool SpatiocyteWorld::has_species_exact(const Species &sp) const
// {
//     return (*space_).has_species_exact(sp);
// }

Integer SpatiocyteWorld::num_molecules(const Species& sp) const
{
    return (*space_).num_molecules(sp);
}

Integer SpatiocyteWorld::num_molecules_exact(const Species& sp) const
{
    return (*space_).num_molecules_exact(sp);
}

Integer SpatiocyteWorld::num_particles(const Species& sp) const
{
    return (*space_).num_particles(sp);
}

Integer SpatiocyteWorld::num_particles_exact(const Species& sp) const
{
    return (*space_).num_particles_exact(sp);
}

Integer SpatiocyteWorld::num_particles() const
{
    return (*space_).num_particles();
}

Integer SpatiocyteWorld::num_voxels() const
{
    return (*space_).num_voxels();
}

Integer SpatiocyteWorld::num_voxels(const Species& sp) const
{
    return (*space_).num_voxels(sp);
}

Integer SpatiocyteWorld::num_voxels_exact(const Species& sp) const
{
    return (*space_).num_voxels_exact(sp);
}

bool SpatiocyteWorld::has_particle(const ParticleID& pid) const
{
    return (*space_).has_particle(pid);
}

bool SpatiocyteWorld::has_voxel(const ParticleID& pid) const
{
    return (*space_).has_voxel(pid);
}

std::vector<std::pair<ParticleID, Particle> >
SpatiocyteWorld::list_particles() const
{
    return (*space_).list_particles();
}

std::vector<std::pair<ParticleID, Particle> >
SpatiocyteWorld::list_particles(const Species& sp) const
{
    return (*space_).list_particles(sp);
}

std::vector<std::pair<ParticleID, Particle> >
SpatiocyteWorld::list_particles_exact(const Species& sp) const
{
    return (*space_).list_particles_exact(sp);
}

std::vector<std::pair<ParticleID, Particle> >
SpatiocyteWorld::list_structure_particles() const
{
    const std::vector<Species> structure_species(list_structure_species());

    typedef std::vector<std::vector<std::pair<ParticleID, Particle> > > tmp_type;
    tmp_type tmp_vector(structure_species.size());
    Integer num_elements;

    for (std::vector<Species>::const_iterator itr(structure_species.begin());
            itr != structure_species.end(); ++itr)
    {
        std::vector<std::pair<ParticleID, Particle> > tmp(list_particles(*itr));
        tmp_vector.push_back(tmp);
        num_elements += tmp.size();
    }

    std::vector<std::pair<ParticleID, Particle> > retval;
    retval.reserve(num_elements);
    for (tmp_type::const_iterator itr(tmp_vector.begin());
            itr != tmp_vector.end(); ++itr)
    {
        retval.insert(retval.end(), (*itr).begin(), (*itr).end());
    }

    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
SpatiocyteWorld::list_non_structure_particles() const
{
    const std::vector<Species> non_structure_species(list_non_structure_species());

    typedef std::vector<std::vector<std::pair<ParticleID, Particle> > > tmp_type;
    tmp_type tmp_vector(non_structure_species.size());
    Integer num_elements;

    for (std::vector<Species>::const_iterator itr(non_structure_species.begin());
            itr != non_structure_species.end(); ++itr)
    {
        std::vector<std::pair<ParticleID, Particle> > tmp(list_particles(*itr));
        tmp_vector.push_back(tmp);
        num_elements += tmp.size();
    }

    std::vector<std::pair<ParticleID, Particle> > retval;
    retval.reserve(num_elements);
    for (tmp_type::const_iterator itr(tmp_vector.begin());
            itr != tmp_vector.end(); ++itr)
    {
        retval.insert(retval.end(), (*itr).begin(), (*itr).end());
    }

    return retval;
}

// bool SpatiocyteWorld::update_particle(const ParticleID& pid, const Particle& p)
// {
//     return (*space_).update_particle(pid, p);
// }

std::vector<Species> SpatiocyteWorld::list_species() const
{
    return (*space_).list_species();
}

std::vector<Species> SpatiocyteWorld::list_non_structure_species() const
{
    const std::vector<Species> species(list_species());
    std::vector<Species> retval;
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        if (!find_molecular_type(*itr)->is_structure())
            retval.push_back(*itr);
    }
    return retval;
}

std::vector<Species> SpatiocyteWorld::list_structure_species() const
{
    const std::vector<Species> species(list_species());
    std::vector<Species> retval;
    for (std::vector<Species>::const_iterator itr(species.begin());
            itr != species.end(); ++itr)
    {
        if (find_molecular_type(*itr)->is_structure())
            retval.push_back(*itr);
    }
    return retval;
}

std::vector<std::pair<ParticleID, Voxel> > SpatiocyteWorld::list_voxels() const
{
    return (*space_).list_voxels();
}

std::vector<std::pair<ParticleID, Voxel> >
    SpatiocyteWorld::list_voxels(const Species& sp) const
{
    return (*space_).list_voxels(sp);
}

std::vector<std::pair<ParticleID, Voxel> >
    SpatiocyteWorld::list_voxels_exact(const Species& sp) const
{
    return (*space_).list_voxels_exact(sp);
}

MolecularTypeBase* SpatiocyteWorld::find_molecular_type(const Species& species)
{
    return (*space_).find_molecular_type(species);
}

const MolecularTypeBase* SpatiocyteWorld::find_molecular_type(const Species& species) const
{
    return (*space_).find_molecular_type(species);
}

MolecularTypeBase* SpatiocyteWorld::get_molecular_type_private(
        const private_coordinate_type& coord)
{
    return (*space_).get_molecular_type(coord);
}

std::pair<std::pair<ParticleID, Voxel>, bool>
SpatiocyteWorld::new_voxel(const Voxel& v)
{
    const private_coordinate_type private_coord(coord2private(v.coordinate()));
    return new_voxel_private(
        Voxel(v.species(), private_coord, v.radius(), v.D(), v.loc()));
}

std::pair<std::pair<ParticleID, Voxel>, bool>
SpatiocyteWorld::new_voxel(const Species& sp, const coordinate_type& coord)
{
    const private_coordinate_type private_coord(coord2private(coord));
    const molecule_info_type minfo(get_molecule_info(sp));
    return new_voxel_private(
        Voxel(sp, private_coord, minfo.radius, minfo.D, minfo.loc));
}

std::pair<std::pair<ParticleID, Voxel>, bool>
SpatiocyteWorld::new_voxel_private(
        const Species& sp, const private_coordinate_type& coord)
{
    const molecule_info_type minfo(get_molecule_info(sp));
    return new_voxel_private(
        Voxel(sp, coord, minfo.radius, minfo.D, minfo.loc));
}

std::pair<std::pair<ParticleID, Voxel>, bool>
SpatiocyteWorld::new_voxel_private(const Voxel& v)
{
    ParticleID pid(sidgen_());
    const bool is_succeeded((*space_).update_voxel_private(pid, v));
    const coordinate_type coord(private2coord(v.coordinate()));
    return std::make_pair(
        std::make_pair(pid, Voxel(v.species(), coord, v.radius(), v.D(), v.loc())),
        is_succeeded);
}

std::pair<std::pair<ParticleID, Voxel>, bool>
SpatiocyteWorld::new_voxel_structure(const Species& sp, const coordinate_type& coord)
{
    const private_coordinate_type private_coord(coord2private(coord));
    const molecule_info_type minfo(get_molecule_info(sp));
    return new_voxel_structure_private(
        Voxel(sp, private_coord, minfo.radius, minfo.D, minfo.loc));
}

std::pair<std::pair<ParticleID, Voxel>, bool>
SpatiocyteWorld::new_voxel_structure_private(const Voxel& v)
{
    const bool is_succeeded((*space_).update_voxel_private(ParticleID(), v));
    const coordinate_type coord(private2coord(v.coordinate()));
    return std::make_pair(std::make_pair(ParticleID(),
                Voxel(v.species(), coord, v.radius(), v.D(), v.loc())),
        is_succeeded);
}

bool SpatiocyteWorld::add_molecules(const Species& sp, const Integer& num)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));

    Integer count(0);
    while (count < num)
    {
        const coordinate_type coord(rng()->uniform_int(0, (*space_).size() - 1));
        const Voxel v(sp, coord2private(coord), info.radius, info.D, info.loc);

        if ((*space_).on_structure(v))
        {
            continue;
        }
        else if (new_voxel_private(v).second)
        {
            ++count;
        }
    }
    return true;
}

bool SpatiocyteWorld::add_molecules(
    const Species& sp, const Integer& num, const boost::shared_ptr<const Shape> shape)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));

    Integer count(0);
    while (count < num)
    {
        const Real3 pos(shape->draw_position(rng_));
        const Voxel v(sp, (*space_).position2private(pos), info.radius, info.D, info.loc);

        if ((*space_).on_structure(v))
        {
            continue;
        }
        else if (new_voxel_private(v).second)
        {
            ++count;
        }
    }
    return true;
}

Integer SpatiocyteWorld::add_structure(
    const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));
    (*space_).make_structure_type(sp, shape->dimension(), info.loc);

    switch (shape->dimension())
    {
    case Shape::THREE:
        return add_structure3(sp, shape);
    case Shape::TWO:
        return add_structure2(sp, shape);
    case Shape::ONE:
    case Shape::UNDEF:
        break;
    }

    throw NotSupported("The dimension of a shape must be two or three.");
}

Integer SpatiocyteWorld::add_structure3(const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    // Real3 l, u;
    // shape->bounding_box(edge_lengths(), l, u);
    // const Real sigma(voxel_radius() * 2);
    // const unsigned int ndim(3);
    // for (unsigned int i(0); i != ndim; ++i)
    // {
    //     l[i] = std::max(0.0, l[i] - sigma);
    //     u[i] = std::min(edge_lengths()[i], u[i] + sigma);
    // }
    // const Integer3 lower(position2global(l)), upper(position2global(u));

    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));
    Integer count(0);
    for (Integer col(0); col < col_size(); ++col)
    {
        for (Integer row(0); row < row_size(); ++row)
        {
            for (Integer layer(0); layer < layer_size(); ++layer)
            {
    // for (Integer col(lower[0]); col < upper[0]; ++col)
    // {
    //     for (Integer row(lower[1]); row < upper[1]; ++row)
    //     {
    //         for (Integer layer(lower[2]); layer < upper[2]; ++layer)
    //         {
                const Integer3 g(col, row, layer);
                const Real L(shape->is_inside(global2position(g)));
                if (L > 0)
                {
                    continue;
                }

                const Voxel v(sp, (*space_).global2private_coord(g),
                    info.radius, info.D, info.loc);
                if (new_voxel_structure_private(v).second)
                {
                    ++count;
                }
            }
        }
    }
    return count;
}

Integer SpatiocyteWorld::add_structure2(const Species& sp, const boost::shared_ptr<const Shape> shape)
{
    // std::ofstream fout("shape.csv");
    // fout << "# " << sp.serial() << std::endl;
    // const unsigned int n(50);
    // const Real3 L(edge_lengths() / static_cast<Real>(n));
    // for (unsigned int i(0); i < n * n * n; ++i)
    // {
    //     const unsigned int x(i % n);
    //     const unsigned int y((i / n) % n);
    //     const unsigned int z(i / (n * n));
    //     if (shape->test_AABB(Real3(x * L[0], y * L[1], z * L[2]),
    //         Real3((x + 1) * L[0], (y + 1) * L[1], (z + 1) * L[2])))
    //     {
    //         fout << x << "," << y << "," << z << std::endl;
    //     }
    // }
    // fout.close();

    // Real3 l, u;
    // shape->bounding_box(edge_lengths(), l, u);
    // const Real sigma(voxel_radius() * 2);
    // const unsigned int ndim(3);
    // for (unsigned int i(0); i != ndim; ++i)
    // {
    //     l[i] = std::max(0.0, l[i] - sigma);
    //     u[i] = std::min(edge_lengths()[i], u[i] + sigma);
    // }
    // const Integer3 lower(position2global(l)), upper(position2global(u));

    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));
    Integer count(0);
    for (Integer col(0); col < col_size(); ++col)
    {
        for (Integer row(0); row < row_size(); ++row)
        {
            for (Integer layer(0); layer < layer_size(); ++layer)
            {
    // for (Integer col(lower[0]); col < upper[0]; ++col)
    // {
    //     for (Integer row(lower[1]); row < upper[1]; ++row)
    //     {
    //         for (Integer layer(lower[2]); layer < upper[2]; ++layer)
    //         {
                const Integer3 g(col, row, layer);
                if (!is_surface_voxel(g, shape))
                {
                    continue;
                }

                const Voxel v(sp, (*space_).global2private_coord(g),
                    info.radius, info.D, info.loc);
                if (new_voxel_structure_private(v).second)
                {
                    ++count;
                }
            }
        }
    }
    return count;
}

bool SpatiocyteWorld::is_surface_voxel(
    const Integer3& g, const boost::shared_ptr<const Shape> shape) const
{
    const Real L(shape->is_inside(global2position(g)));
    if (L > 0 || L < -2 * voxel_radius())
    {
        return false;
    }

    const SpatiocyteWorld::private_coordinate_type
        private_coord((*space_).global2private_coord(g));
    for (Integer i(0); i < 12; ++i)
    {
        if (shape->is_inside(global2position((*space_).private_coord2global(
            (*space_).get_neighbor_private(private_coord, i)))) > 0)
        {
            return true;
        }
    }
    return false;
}

// TODO
Integer SpatiocyteWorld::add_neighbors(const Species& sp,
    const SpatiocyteWorld::private_coordinate_type center)
{
    Integer count(0);
    const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));
    for (Integer i(0); i < 12; ++i)
    {
        const private_coordinate_type n((*space_).get_neighbor_private(center, i));
        if (new_voxel_private(Voxel(sp, n, info.radius, info.D, info.loc)).second)
        {
            ++count;
        }
        else
        {
            throw "Error in add_neighbors()";
        }
    }
    return count;

    // Integer count(0);
    // const SpatiocyteWorld::molecule_info_type info(get_molecule_info(sp));
    // std::vector<SpatiocyteWorld::private_coordinate_type> neighbors(
    //         (*space_).get_neighbors(center));
    // for (std::vector<SpatiocyteWorld::private_coordinate_type>::iterator itr(
    //             neighbors.begin()); itr != neighbors.end(); itr++)
    // {
    //     if (new_voxel_private(Voxel(sp, *itr, info.radius, info.D, info.loc)).second)
    //         ++count;
    //     else
    //         throw "Error in add_neighbors()";
    // }
    // return count;
}
// TODO

void SpatiocyteWorld::remove_molecules(const Species& sp, const Integer& num)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }

    MolecularTypeBase* mtype(find_molecular_type(sp));
    if (!mtype->with_voxels())
    {
        throw NotSupported(
            "remove_molecuels for MolecularType with no voxel is not supported now");
    }
    else if (mtype->size() < num)
    {
        throw std::invalid_argument(
            "The number of molecules cannot be negative.");
    }

    Integer count(0);
    while (count < num)
    {
        const Integer idx(rng_->uniform_int(0, mtype->size() - 1));
        if (remove_voxel_private(mtype->at(idx).first))
        {
            ++count;
        }
    }
}

bool SpatiocyteWorld::remove_voxel_private(const private_coordinate_type coord)
{
    return (*space_).remove_voxel_private(coord);
}

bool SpatiocyteWorld::move(coordinate_type from, coordinate_type to)
{
    return (*space_).move(from, to);
}

bool SpatiocyteWorld::move_private(const private_coordinate_type& src,
        const private_coordinate_type& dest, const std::size_t candidate)
{
    return (*space_).move_private(src, dest, candidate);
}

bool SpatiocyteWorld::can_move(const private_coordinate_type& src,
        const private_coordinate_type& dest) const
{
    return (*space_).can_move(src, dest);
}

std::pair<SpatiocyteWorld::private_coordinate_type, bool>
SpatiocyteWorld::move_to_neighbor(
    MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
    particle_info_type& info, const Integer nrand)
{
    return (*space_).move_to_neighbor(from_mt, loc, info, nrand);
}

std::pair<SpatiocyteWorld::private_coordinate_type, bool>
SpatiocyteWorld::check_neighbor_private(
    const private_coordinate_type coord, const std::string& loc)
{
    std::vector<private_coordinate_type> tmp;
    tmp.reserve(12);
    for (unsigned int rnd(0); rnd < 12; ++rnd)
    {
        const private_coordinate_type
            neighbor((*space_).get_neighbor_private(coord, rnd));
        const MolecularTypeBase* mt(get_molecular_type_private(neighbor));
        const std::string
            serial(mt->is_vacant() ? "" : mt->species().serial());
        if (serial == loc)
        {
            tmp.push_back(neighbor);
        }
    }

    if (tmp.size() == 0)
    {
        return std::make_pair(coord, false);
    }

    return std::make_pair(
        tmp[rng()->uniform_int(0, tmp.size() - 1)], true);

    // const Integer rnd(rng()->uniform_int(0, 11));
    // const private_coordinate_type neighbor((*space_).get_neighbor_private(coord, rnd));
    // bool flg = get_molecular_type_private(neighbor)->is_vacant(); //XXX: loc
    // return std::make_pair(neighbor, flg);
}

} // spatiocyte

} // ecell4
