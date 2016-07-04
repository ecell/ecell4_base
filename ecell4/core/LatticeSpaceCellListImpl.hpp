#ifndef __ECELL4_LATTICE_SPACE_CELL_LIST_IMPL_HPP
#define __ECELL4_LATTICE_SPACE_CELL_LIST_IMPL_HPP

#include "Context.hpp"
#include "MolecularType.hpp"
#include "VacantType.hpp"
// #include <cmath>
#include <sstream>

#include "comparators.hpp"
#include "LatticeSpace.hpp"

namespace ecell4
{

inline unsigned int ceilint(const unsigned int x, const unsigned int y)
{
    return x / y + (x % y != 0);
}

class LatticeSpaceCellListImpl
    : public LatticeSpaceBase
{
public:

    typedef LatticeSpaceBase base_type;

    typedef base_type::particle_info_type particle_info_type;
    typedef base_type::private_coordinate_type private_coordinate_type;
    typedef base_type::private_coordinate_type coordinate_type;

    typedef std::map<Species, boost::shared_ptr<MolecularType> > spmap;
    typedef std::vector<std::pair<MolecularTypeBase*, private_coordinate_type> >
        cell_type;
    typedef std::vector<cell_type> matrix_type;
    typedef std::map<Species, boost::shared_ptr<const Shape> > structure_container_type;

public:

    LatticeSpaceCellListImpl(
        const Real3& edge_lengths, const Real& voxel_radius,
        const Integer3& matrix_sizes, const bool is_periodic = true)
        : base_type(edge_lengths, voxel_radius, is_periodic), is_periodic_(is_periodic),
        matrix_sizes_(matrix_sizes),
        matrix_(matrix_sizes_[0] * matrix_sizes_[1] * matrix_sizes_[2])
    {
        cell_sizes_[0] = ceilint(col_size_, matrix_sizes_[0]);
        cell_sizes_[1] = ceilint(row_size_, matrix_sizes_[1]);
        cell_sizes_[2] = ceilint(layer_size_, matrix_sizes_[2]);

        vacant_ = &(VacantType::getInstance());
        std::stringstream ss;
        ss << voxel_radius_;
        border_ = new MolecularType(Species("Border", ss.str(), "0"));
        periodic_ = new MolecularType(Species("Periodic", ss.str(), "0"));
    }

    virtual ~LatticeSpaceCellListImpl()
    {
        delete border_;
        delete periodic_;
    }

    /**
     */

    inline matrix_type::size_type coord2index(const private_coordinate_type& coord) const
    {
        return global2index(private_coord2global(coord));
    }

    inline matrix_type::size_type global2index(const Integer3& g) const
    {
        return (g.col / cell_sizes_[0]) + matrix_sizes_[0] * ((g.row / cell_sizes_[1]) + matrix_sizes_[1] * (g.layer / cell_sizes_[2]));
    }

    cell_type::iterator find_from_cell(
        const private_coordinate_type& coord, cell_type& cell)
    {
        return std::find_if(cell.begin(), cell.end(),
            utils::pair_second_element_unary_predicator<
                MolecularTypeBase*, private_coordinate_type>(coord));
        // cell_type::iterator i(cell.begin());
        // for (; i != cell.end(); ++i)
        // {
        //     if ((*i).second == coord)
        //     {
        //         return i;
        //     }
        // }
        // return i;
    }

    cell_type::const_iterator find_from_cell(
        const private_coordinate_type& coord, const cell_type& cell) const
    {
        return std::find_if(cell.begin(), cell.end(),
            utils::pair_second_element_unary_predicator<
                MolecularTypeBase*, private_coordinate_type>(coord));
        // cell_type::const_iterator i(cell.begin());
        // for (; i != cell.end(); ++i)
        // {
        //     if ((*i).second == coord)
        //     {
        //         return i;
        //     }
        // }
        // return i;
    }

    void update_matrix(const private_coordinate_type& coord, MolecularTypeBase* mt)
    {
        cell_type& cell(matrix_[coord2index(coord)]);
        cell_type::iterator i(find_from_cell(coord, cell));

        if (i != cell.end())
        {
            if (mt->is_vacant())
            {
                cell.erase(i);
            }
            else
            {
                (*i).first = mt;
            }
        }
        else if (!mt->is_vacant())
        {
            cell.push_back(std::make_pair(mt, coord));
        }
        else
        {
            throw NotFound("1");
        }
    }

    void update_matrix(const private_coordinate_type& from_coord,
        const private_coordinate_type& to_coord,
        MolecularTypeBase* mt)
    {
        const matrix_type::size_type from_idx(coord2index(from_coord)),
            to_idx(coord2index(to_coord));
        if (from_idx == to_idx)
        {
            cell_type& cell(matrix_[from_idx]);
            cell_type::iterator i(find_from_cell(from_coord, cell));
            if (i == cell.end())
            {
                throw NotFound("2");
            }
            (*i).first = mt;
            (*i).second = to_coord;
        }
        else
        {
            cell_type& cell(matrix_[from_idx]);
            cell_type::iterator i(find_from_cell(from_coord, cell));
            if (i == cell.end())
            {
                throw NotFound("3");
            }
            cell.erase(i);
            matrix_[to_idx].push_back(std::make_pair(mt, to_coord));
        }
    }

    void dump_matrix()
    {
        std::cout << "=====================" << std::endl;
        for (matrix_type::size_type i(0); i != matrix_.size(); ++i)
        {
            cell_type& c(matrix_[i]);
            if (c.size() == 0)
            {
                continue;
            }

            std::cout << i << " : ";
            for (cell_type::const_iterator j(c.begin()); j != c.end(); ++j)
            {
                std::cout << (*j).second << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "=====================" << std::endl;
    }

    /**
     */

    virtual std::vector<Species> list_species() const
    {
        std::vector<Species> keys;
        keys.reserve(spmap_.size());
        for (spmap::const_iterator i(spmap_.begin());
            i != spmap_.end(); ++i)
        {
            keys.push_back((*i).first);
        }
        return keys;
    }

    Integer count_voxels(
        const boost::shared_ptr<MolecularType>& mt) const
    {
        Integer count(0);
        utils::pair_first_element_unary_predicator<
            MolecularTypeBase*, private_coordinate_type> pred(mt.get());

        for (matrix_type::const_iterator i(matrix_.begin());
            i != matrix_.end(); ++i)
        {
            count += static_cast<Integer>(
                std::count_if((*i).begin(), (*i).end(), pred));
        }
        return count;
    }

    virtual Integer num_voxels_exact(const Species& sp) const
    {
        spmap::const_iterator itr(spmap_.find(sp));
        if (itr == spmap_.end())
        {
            return 0;
        }
        const boost::shared_ptr<MolecularType>& mt((*itr).second);
        return mt->size();
    }

    virtual Integer num_voxels(const Species& sp) const
    {
        Integer count(0);
        SpeciesExpressionMatcher sexp(sp);
        for (spmap::const_iterator itr(spmap_.begin());
                itr != spmap_.end(); ++itr)
        {
            if (sexp.match((*itr).first))
            {
                const boost::shared_ptr<MolecularType>& mt((*itr).second);
                count += mt->size();
            }
        }
        return count;
    }

    virtual Integer num_voxels() const
    {
        Integer retval(0);
        for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
        {
            retval += (*itr).second->size();
        }
        return retval;
    }

    virtual bool has_voxel(const ParticleID& pid) const
    {
        for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
        {
            const boost::shared_ptr<MolecularType>& mt((*itr).second);
            if (mt->is_vacant())
            {
                return false;
            }
            for (MolecularType::const_iterator vitr(mt->begin());
                vitr != mt->end(); ++vitr)
            {
                if ((*vitr).second == pid)
                {
                    return true;
                }
            }
        }
        return false;
    }

    virtual std::vector<std::pair<ParticleID, Voxel> > list_voxels() const
    {
        std::vector<std::pair<ParticleID, Voxel> > retval;

        for (spmap::const_iterator itr(spmap_.begin()); itr != spmap_.end(); ++itr)
        {
            const boost::shared_ptr<MolecularType>& mt((*itr).second);
            const std::string loc((mt->location()->is_vacant())
                ? "" : mt->location()->species().serial());
            const Species& sp(mt->species());
            for (MolecularType::const_iterator itr(mt->begin());
                itr != mt->end(); ++itr)
            {
                retval.push_back(std::make_pair(
                    (*itr).second,
                    Voxel(sp, private2coord((*itr).first), mt->radius(), mt->D(), loc)));
            }
        }
        return retval;
    }

    virtual std::vector<std::pair<ParticleID, Voxel> >
        list_voxels(const Species& sp) const
    {
        SpeciesExpressionMatcher sexp(sp);
        std::vector<std::pair<ParticleID, Voxel> > retval;
        for (spmap::const_iterator itr(spmap_.begin());
                itr != spmap_.end(); ++itr)
        {
            if (!sexp.match((*itr).first))
            {
                continue;
            }

            const boost::shared_ptr<MolecularType>& mt((*itr).second);
            const std::string loc((mt->location()->is_vacant())
                ? "" : mt->location()->species().serial());
            for (MolecularType::const_iterator i(mt->begin());
                i != mt->end(); ++i)
            {
                retval.push_back(std::make_pair(
                    (*i).second,
                    Voxel(sp, private2coord((*i).first), mt->radius(), mt->D(), loc)));
            }
        }
        return retval;
    }

    virtual std::vector<std::pair<ParticleID, Voxel> >
        list_voxels_exact(const Species& sp) const
    {
        std::vector<std::pair<ParticleID, Voxel> > retval;
        spmap::const_iterator itr(spmap_.find(sp));
        if (itr == spmap_.end())
        {
            return retval;
        }

        const boost::shared_ptr<MolecularType>& mt((*itr).second);
        const std::string loc((mt->location()->is_vacant())
            ? "" : mt->location()->species().serial());
        for (MolecularType::const_iterator i(mt->begin());
            i != mt->end(); ++i)
        {
            retval.push_back(std::make_pair(
                (*i).second,
                Voxel(sp, private2coord((*i).first), mt->radius(), mt->D(), loc)));
        }
        return retval;
    }

    /**
     * Change the Species at v.coordinate() to v.species.
     * The ParticleID must be kept after this update.
     */
    virtual void update_voxel_private(const Voxel& v)
    {
        const private_coordinate_type coord(v.coordinate());
        MolecularTypeBase* src_mt(get_molecular_type(coord));
        MolecularTypeBase* new_mt(get_molecular_type(v));

        if (src_mt->with_voxels() != new_mt->with_voxels())
        {
            throw NotSupported("ParticleID is needed/lost.");
        }

        new_mt->add_voxel_without_checking(src_mt->pop(coord));
        update_matrix(coord, new_mt);

        // const private_coordinate_type coord(v.coordinate());
        // MolecularTypeBase* src_mt(get_molecular_type(coord));
        // if (src_mt->is_vacant())
        // {
        //     return false;
        // }

        // MolecularTypeBase* new_mt(get_molecular_type(v));
        // if (new_mt->is_vacant())
        // {
        //     // ???
        //     return false;
        // }

        // new_mt->add_voxel_without_checking(src_mt->pop(coord));
        // update_matrix(coord, new_mt);
        // return true;
    }

    virtual bool update_voxel_private(const ParticleID& pid, const Voxel& v);

    virtual std::pair<ParticleID, Voxel> get_voxel(const ParticleID& pid) const
    {
        const std::pair<const MolecularTypeBase*, private_coordinate_type>
            target(__get_coordinate(pid));
        if (target.second == -1)
        {
            throw NotFound("voxel not found.");
        }

        const coordinate_type coord(private2coord(target.second));
        const MolecularTypeBase* mt(target.first);
        const std::string loc((mt->location()->is_vacant())
            ? "" : mt->location()->species().serial());
        return std::make_pair(
            pid, Voxel(mt->species(), coord, mt->radius(), mt->D(), loc));
    }

    virtual std::pair<ParticleID, Voxel> get_voxel(const coordinate_type& coord) const
    {
        const private_coordinate_type private_coord(coord2private(coord));
        const MolecularTypeBase* mt(get_molecular_type(private_coord));
        const std::string loc((mt->location()->is_vacant())
            ? "" : mt->location()->species().serial());
        if (mt->with_voxels())
        {
            return std::make_pair(mt->find_particle_id(private_coord),
                Voxel(mt->species(), coord, mt->radius(), mt->D(), loc));
        }
        else
        {
            return std::make_pair(ParticleID(),
                Voxel(mt->species(), coord, mt->radius(), mt->D(), loc));
        }
    }

    virtual std::pair<ParticleID, Voxel> get_voxel_private(const coordinate_type& private_coord) const
    {
        const MolecularTypeBase* mt(get_molecular_type(private_coord));
        const std::string loc((mt->location()->is_vacant())
            ? "" : mt->location()->species().serial());
        if (mt->with_voxels())
        {
            return std::make_pair(mt->find_particle_id(private_coord),
                Voxel(mt->species(), private_coord, mt->radius(), mt->D(), loc));
        }
        else
        {
            return std::make_pair(ParticleID(),
                Voxel(mt->species(), private_coord, mt->radius(), mt->D(), loc));
        }
    }

    virtual bool remove_voxel(const ParticleID& pid)
    {
        std::pair<MolecularTypeBase*, private_coordinate_type>
            target(__get_coordinate(pid));
        if (target.second != -1)
        {
            MolecularTypeBase* mt(target.first);
            const private_coordinate_type coord(target.second);
            if (!mt->remove_voxel_if_exists(coord))
            {
                return false;
            }

            mt->location()->add_voxel_without_checking(
                particle_info_type(coord, ParticleID()));
            update_matrix(coord, mt->location());
            return true;
        }
        return false;
    }

    virtual bool remove_voxel_private(const private_coordinate_type& coord)
    {
        MolecularTypeBase* mt(get_molecular_type(coord));
        if (mt->is_vacant())
        {
            return false;
        }

        if (mt->remove_voxel_if_exists(coord))
        {
            // ???
            update_matrix(coord, vacant_);
            return true;
        }
        return true;
    }

    virtual bool move(const coordinate_type& from, const coordinate_type& to)
    {
        const private_coordinate_type src(coord2private(from));
        const private_coordinate_type dest(coord2private(to));
        return move_private(src, dest);
    }

    virtual bool move_private(const private_coordinate_type& src,
            const private_coordinate_type& dest, const std::size_t candidate=0)
    {
        private_coordinate_type tmp_dest(dest);
        if (src == tmp_dest)
        {
            return false;
        }

        MolecularTypeBase* src_mt(get_molecular_type(src));
        if (src_mt->is_vacant())
        {
            return true;
        }

        MolecularTypeBase* dest_mt(get_molecular_type(tmp_dest));
        if (dest_mt == border_)
        {
            return false;
        }
        else if (dest_mt == periodic_)
        {
            tmp_dest = periodic_transpose_private(tmp_dest);
            dest_mt = get_molecular_type(tmp_dest);
        }

        if (dest_mt != src_mt->location())
        {
            return false;
        }

        src_mt->replace_voxel(src, tmp_dest);
        dest_mt->replace_voxel(tmp_dest, src);
        if (!dest_mt->is_vacant())
        {
            update_matrix(src, dest_mt);
            update_matrix(tmp_dest, src_mt);
        }
        else
        {
            update_matrix(src, tmp_dest, src_mt);
        }
        return true;
    }

    virtual const Particle particle_at(const coordinate_type& coord) const
    {
        const private_coordinate_type private_coord(coord2private(coord));
        const MolecularTypeBase* mt(get_molecular_type(private_coord));
        return Particle(
            mt->species(), coordinate2position(coord), mt->radius(), mt->D());
    }

    virtual MolecularTypeBase* find_molecular_type(const Species& sp)
    {
        spmap::iterator itr(spmap_.find(sp));
        if (itr == spmap_.end())
        {
            throw NotFound("MolecularType not found.");
        }
        return (*itr).second.get(); //XXX: Raw pointer was thrown.
    }

    virtual const MolecularTypeBase* find_molecular_type(const Species& sp) const
    {
        spmap::const_iterator itr(spmap_.find(sp));
        if (itr == spmap_.end())
        {
            throw NotFound("MolecularType not found.");
        }
        return (*itr).second.get(); //XXX: Raw pointer was thrown.
    }

    virtual MolecularTypeBase* get_molecular_type(
        const private_coordinate_type& coord);
    const MolecularTypeBase* get_molecular_type(
        const private_coordinate_type& coord) const;
    MolecularTypeBase* get_molecular_type(const Voxel& v);

    virtual bool on_structure(const Voxel& v)
    {
        return (get_molecular_type(v.coordinate())
            != get_molecular_type(v)->location()); //XXX: == ???
    }

    private_coordinate_type get_neighbor_private_boundary(
        const private_coordinate_type& coord, const Integer& nrand) const
    {
        private_coordinate_type const dest = get_neighbor_private(coord, nrand);
        return (!is_periodic_ || is_inside(dest) ? dest : periodic_transpose_private(dest));
    }

    virtual void add_structure(const Species& sp,
        const boost::shared_ptr<const Shape>& s, const std::string loc);
    virtual const boost::shared_ptr<const Shape>& get_structure(const Species& sp) const;
    virtual const Shape::dimension_kind get_structure_dimension(const Species& sp) const;

    virtual std::pair<private_coordinate_type, bool> move_to_neighbor(
        MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
        particle_info_type& info, const Integer nrand);

    virtual Integer num_molecules(const Species& sp) const;

    /**
     */

    virtual Real get_value(const Species& sp) const
    {
        return static_cast<Real>(num_molecules(sp));
    }

    virtual Real get_value_exact(const Species& sp) const
    {
        return static_cast<Real>(num_molecules_exact(sp));
    }

    bool has_species(const Species& sp) const
    {
        return spmap_.find(sp) != spmap_.end();
    }

protected:

    std::pair<spmap::iterator, bool> __get_molecular_type(const Voxel& v);
    std::pair<MolecularTypeBase*, private_coordinate_type>
        __get_coordinate(const ParticleID& pid);
    std::pair<const MolecularTypeBase*, private_coordinate_type>
        __get_coordinate(const ParticleID& pid) const;

    bool make_structure_type(
        const Species& sp, Shape::dimension_kind dimension, const std::string loc);

protected:

    bool is_periodic_;

    spmap spmap_;

    MolecularTypeBase* vacant_;
    MolecularTypeBase* border_;
    MolecularTypeBase* periodic_;

    Integer3 matrix_sizes_, cell_sizes_;
    matrix_type matrix_;
    structure_container_type structures_;
};

} // ecell4

#endif /* __ECELL4_LATTICE_SPACE_CELL_LIST_IMPL_HPP */
