#include <numeric>
#include "SubvolumeSpace.hpp"
#include "Context.hpp"

namespace ecell4
{

Integer SubvolumeSpaceVectorImpl::num_molecules(const Species& sp) const
{
    SpeciesExpressionMatcher sexp(sp);
    Integer retval(0);
    for (matrix_type::const_iterator i(matrix_.begin());
        i != matrix_.end(); ++i)
    {
        if (sexp.match((*i).first))
        {
            do
            {
                // retval += std::accumulate((*i).second.begin(), (*i).second.end(), 0);
                retval += (*i).second->num_molecules();
            } while (sexp.next());
        }
    }
    return retval;
}

Integer SubvolumeSpaceVectorImpl::num_molecules_exact(const Species& sp) const
{
    matrix_type::const_iterator i(matrix_.find(sp));
    if (i == matrix_.end())
    {
        return 0;
    }
    // return std::accumulate((*i).second.begin(), (*i).second.end(), 0);
    return (*i).second->num_molecules();
}

Integer SubvolumeSpaceVectorImpl::num_molecules(
    const Species& sp, const coordinate_type& c) const
{
    SpeciesExpressionMatcher sexp(sp);
    Integer retval(0);
    for (matrix_type::const_iterator i(matrix_.begin());
        i != matrix_.end(); ++i)
    {
        if (sexp.match((*i).first))
        {
            do
            {
                // retval += (*i).second[c];
                retval += (*i).second->num_molecules(c);
            } while (sexp.next());
        }
    }
    return retval;
}

Integer SubvolumeSpaceVectorImpl::num_molecules_exact(
    const Species& sp, const coordinate_type& c) const
{
    matrix_type::const_iterator i(matrix_.find(sp));
    if (i == matrix_.end())
    {
        return 0;
    }
    // return (*i).second[c];
    return (*i).second->num_molecules(c);
}

const boost::shared_ptr<SubvolumeSpaceVectorImpl::PoolBase>&
SubvolumeSpaceVectorImpl::get_pool(const Species& sp) const
{
    matrix_type::const_iterator i(matrix_.find(sp));
    if (i == matrix_.end())
    {
        std::ostringstream message;
        message << "Speices [" << sp.serial() << "] not found";
        throw NotFound(message.str());
        // return boost::shared_ptr<PoolBase>();
    }
    return (*i).second;
}

const boost::shared_ptr<SubvolumeSpaceVectorImpl::PoolBase>
SubvolumeSpaceVectorImpl::reserve_pool(
    const Species& sp, const Real D, const Species::serial_type& loc)
{
    matrix_type::const_iterator i(matrix_.find(sp));
    if (i != matrix_.end())
    {
        throw AlreadyExists("Species already exists");
    }
    // matrix_.insert(std::make_pair(sp, std::vector<Integer>(num_subvolumes())));
    const boost::shared_ptr<PoolBase> pool(new Pool(sp, D, loc, num_subvolumes()));
    matrix_.insert(std::make_pair(sp, pool));
    species_.push_back(sp);
    return pool;
}

void SubvolumeSpaceVectorImpl::add_molecules(
    const Species& sp, const Integer& num, const coordinate_type& c)
{
    matrix_type::iterator i(matrix_.find(sp));
    if (i == matrix_.end())
    {
        if (has_structure(sp))
        {
            return;
        }

        // reserve_species(sp, c);
        // i = matrix_.find(sp);
        std::ostringstream message;
        message << "Speices [" << sp.serial() << "] not found";
        throw NotFound(message.str());
    }
    // (*i).second[c] += num;
    (*i).second->add_molecules(num, c);
}

void SubvolumeSpaceVectorImpl::remove_molecules(
    const Species& sp, const Integer& num, const coordinate_type& c)
{
    matrix_type::iterator i(matrix_.find(sp));
    if (i == matrix_.end())
    {
        if (has_structure(sp))
        {
            return;
        }

        std::ostringstream message;
        message << "Speices [" << sp.serial() << "] not found";
        throw NotFound(message.str());
    }

    // if ((*i).second[c] < num)
    if ((*i).second->num_molecules(c) < num)
    {
        std::ostringstream message;
        message << "The number of molecules cannot be negative. [" << sp.serial() << "]";
        throw std::invalid_argument(message.str());
    }

    // (*i).second[c] -= num;
    (*i).second->remove_molecules(num, c);
}

SubvolumeSpaceVectorImpl::coordinate_type SubvolumeSpaceVectorImpl::get_neighbor(
    const coordinate_type& c, const Integer rnd) const
{
    Integer3 g(coord2global(c));

    switch (rnd)
    {
    case 0:
        return global2coord(g.east());
    case 1:
        return global2coord(g.west());
    case 2:
        return global2coord(g.south());
    case 3:
        return global2coord(g.north());
    case 4:
        return global2coord(g.dorsal());
    case 5:
        return global2coord(g.ventral());
    }

    throw IllegalState("the number of neighbors is less than 6.");
}

std::vector<SubvolumeSpaceVectorImpl::coordinate_type>
SubvolumeSpaceVectorImpl::list_coordinates(const Species& sp) const
{
    SpeciesExpressionMatcher sexp(sp);
    std::vector<coordinate_type> retval;
    for (matrix_type::const_iterator i(matrix_.begin());
        i != matrix_.end(); ++i)
    {
        Integer cnt(sexp.count((*i).first));
        if (cnt > 0)
        {
            const std::vector<coordinate_type> coords = (*i).second->list_coordinates();
            for (; cnt > 0; --cnt)
            {
                retval.insert(retval.end(), coords.begin(), coords.end());
            }
            // for (cell_type::size_type j(0); j < (*i).second.size(); ++j)
            // {
            //     if ((*i).second[j] > 0)
            //     {
            //         retval.resize(retval.size() + (*i).second[j] * cnt, j);
            //     }
            // }
        }
    }
    return retval;
}

std::vector<SubvolumeSpaceVectorImpl::coordinate_type>
SubvolumeSpaceVectorImpl::list_coordinates_exact(const Species& sp) const
{
    std::vector<coordinate_type> retval;
    matrix_type::const_iterator i(matrix_.find(sp));
    if (i == matrix_.end())
    {
        return retval;
    }
    return (*i).second->list_coordinates();

    // for (cell_type::size_type j(0); j < (*i).second.size(); ++j)
    // {
    //     if ((*i).second[j] > 0)
    //     {
    //         retval.resize(retval.size() + (*i).second[j], j);
    //     }
    // }
    // return retval;
}

void SubvolumeSpaceVectorImpl::add_structure(
    const Species& sp, const boost::shared_ptr<const Shape>& shape)
{
    structure_matrix_type::const_iterator it(structure_matrix_.find(sp.serial()));
    if (it != structure_matrix_.end())
    {
        std::ostringstream message;
        message << "The given structure [" << sp.serial() << "] is already defined.";
        throw AlreadyExists(message.str());
    }

    switch (shape->dimension())
    {
    case Shape::THREE:
        add_structure3(sp, shape);
        return;
    case Shape::TWO:
        add_structure2(sp, shape);
        return;
    case Shape::ONE:
    case Shape::UNDEF:
        break;
    }

    throw NotSupported("The dimension of a shape must be two or three.");
}

void SubvolumeSpaceVectorImpl::add_structure3(
    const Species& sp, const boost::shared_ptr<const Shape>& shape)
{
    structure_cell_type overlap(num_subvolumes());
    for (structure_cell_type::size_type i(0); i != overlap.size(); ++i)
    {
        if (shape->is_inside(coord2position(i)) > 0)
        {
            overlap[i] = 0;
        }
        else
        {
            overlap[i] = 1;
        }
    }
    // structures_.insert(std::make_pair(sp.serial(), Shape::THREE));
    structure_matrix_.insert(std::make_pair(sp.serial(), overlap));
}

void SubvolumeSpaceVectorImpl::add_structure2(
    const Species& sp, const boost::shared_ptr<const Shape>& shape)
{
    structure_cell_type overlap(num_subvolumes());
    for (structure_cell_type::size_type i(0); i != overlap.size(); ++i)
    {
        if (is_surface_subvolume(i, shape))
        {
            // overlap[i] = 1;
            overlap[i] = unit_area();
        }
        else
        {
            overlap[i] = 0;
        }
    }
    // structures_.insert(std::make_pair(sp.serial(), Shape::TWO));
    structure_matrix_.insert(std::make_pair(sp.serial(), overlap));
}

bool SubvolumeSpaceVectorImpl::is_surface_subvolume(
    const coordinate_type& c, const boost::shared_ptr<const Shape>& shape)
{
    const Real3 lengths(subvolume_edge_lengths());
    const Real3 center(coord2position(c));

    if (shape->is_inside(center) > 0)
    {
        return false;
    }

    for (unsigned int dim(0); dim < 3 * 3 * 3; ++dim)
    {
        const int x(static_cast<int>(dim / 9));
        const int y(static_cast<int>((dim - x * 9) / 3));
        const int z(dim - (x * 3 + y) * 3);

        if ((x == 1 && y == 1 && z == 1)
            || (x != 1 && y != 1 && z != 1))
        {
            continue;
        }

        const Real3 shift(
            (x - 1) * lengths[0], (y - 1) * lengths[1], (z - 1) * lengths[2]);
        const Real3 neighbor = center + shift;
        if (shape->is_inside(neighbor) > 0)
        {
            return true;
        }
    }
    return false;
}

bool SubvolumeSpaceVectorImpl::check_structure(
    const Species::serial_type& serial,
    const SubvolumeSpaceVectorImpl::coordinate_type& coord) const
{
    if (serial == "")
    {
        ; // This is for default structure.
        return true;
    }

    structure_matrix_type::const_iterator i(structure_matrix_.find(serial));
    if (i == structure_matrix_.end())
    {
        return false;
    }
    return ((*i).second[coord] > 0);
}

void SubvolumeSpaceVectorImpl::update_structure(
    const Species::serial_type& serial, const coordinate_type& coord,
    const Real& value)
{
    structure_matrix_type::iterator i(structure_matrix_.find(serial));
    if (i == structure_matrix_.end())
    {
        structure_cell_type overlap(num_subvolumes());
        overlap[coord] = value;
        structure_matrix_.insert(std::make_pair(serial, overlap));
        // structures_.insert(std::make_pair(serial, Shape::THREE));  //XXX: as a default
    }
    else
    {
        (*i).second[coord] = value;
    }
}

std::vector<Species::serial_type> SubvolumeSpaceVectorImpl::list_structures() const
{
    std::vector<Species::serial_type> retval;
    for (structure_matrix_type::const_iterator i(structure_matrix_.begin());
        i != structure_matrix_.end(); ++i)
    {
        retval.push_back((*i).first);
    }
    return retval;
}

Real SubvolumeSpaceVectorImpl::get_volume(const Species& sp) const
{
    structure_matrix_type::const_iterator i(structure_matrix_.find(sp.serial()));
    if (i == structure_matrix_.end())
    {
        return 0.0;
    }
    const Real occupancy = std::accumulate((*i).second.begin(), (*i).second.end(), 0.0);
    return subvolume() * occupancy;
}

const Integer SubvolumeSpaceVectorImpl::num_subvolumes(const Species& sp) const
{
    structure_matrix_type::const_iterator i(structure_matrix_.find(sp.serial()));
    if (i == structure_matrix_.end())
    {
        return 0;
    }

    Integer num(0);
    for (structure_cell_type::const_iterator j((*i).second.begin());
        j != (*i).second.end(); ++j)
    {
        if (*j > 0)
        {
            ++num;
        }
    }
    return num;
    // return std::accumulate((*i).second.begin(), (*i).second.end(), 0);
}

} // ecell4
