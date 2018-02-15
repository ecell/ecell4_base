#ifndef ECELL4_LATTICE_SPACE_CELL_LIST_IMPL_HPP
#define ECELL4_LATTICE_SPACE_CELL_LIST_IMPL_HPP

#include "Context.hpp"
#include "MolecularType.hpp"
#include "VacantType.hpp"
// #include <cmath>
#include <sstream>

#include "comparators.hpp"
#include "HCPLatticeSpace.hpp"

namespace ecell4
{

inline unsigned int ceilint(const unsigned int x, const unsigned int y)
{
    return x / y + (x % y != 0);
}

class LatticeSpaceCellListImpl
    : public HCPLatticeSpace
{
public:

    typedef HCPLatticeSpace base_type;

    typedef base_type::coordinate_id_pair_type coordinate_id_pair_type;
    typedef base_type::coordinate_type coordinate_type;

    typedef std::vector<std::pair<boost::shared_ptr<VoxelPool>, coordinate_type> >
        cell_type;
    typedef std::vector<cell_type> matrix_type;
    typedef std::map<Species, boost::shared_ptr<const Shape> > structure_container_type;

protected:

    typedef utils::get_mapper_mf<
        Species, boost::shared_ptr<VoxelPool> >::type voxel_pool_map_type;
    typedef utils::get_mapper_mf<
        Species, boost::shared_ptr<MoleculePool> >::type molecule_pool_map_type;
    // typedef std::map<
    //     Species, boost::shared_ptr<VoxelPool> > voxel_pool_map_type;
    // typedef std::map<
    //     Species, boost::shared_ptr<MoleculePool> > molecule_pool_map_type;

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

        std::stringstream ss;
        ss << voxel_radius_;
        border_ = boost::shared_ptr<VoxelPool>(
                new MolecularType(Species("Border", ss.str(), "0"), vacant_));
        periodic_ = boost::shared_ptr<VoxelPool>(
                new MolecularType(Species("Periodic", ss.str(), "0"), vacant_));
    }

    virtual ~LatticeSpaceCellListImpl() {}

    /**
     */

    inline matrix_type::size_type coordinate2index(const coordinate_type& coord) const
    {
        return global2index(coordinate2global(coord));
    }

    inline matrix_type::size_type global2index(const Integer3& g) const
    {
        return (g.col / cell_sizes_[0]) + matrix_sizes_[0] * ((g.row / cell_sizes_[1]) + matrix_sizes_[1] * (g.layer / cell_sizes_[2]));
    }

    cell_type::iterator find_from_cell(
        const coordinate_type& coord, cell_type& cell)
    {
        return std::find_if(cell.begin(), cell.end(),
            utils::pair_second_element_unary_predicator<
                boost::shared_ptr<VoxelPool>, coordinate_type>(coord));
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
        const coordinate_type& coord, const cell_type& cell) const
    {
        return std::find_if(cell.begin(), cell.end(),
            utils::pair_second_element_unary_predicator<
                boost::shared_ptr<VoxelPool>, coordinate_type>(coord));
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

    void update_matrix(const coordinate_type& coord, boost::shared_ptr<VoxelPool> vp)
    {
        cell_type& cell(matrix_[coordinate2index(coord)]);
        cell_type::iterator i(find_from_cell(coord, cell));

        if (i != cell.end())
        {
            if (vp->is_vacant())
            {
                cell.erase(i);
            }
            else
            {
                (*i).first = vp;
            }
        }
        else if (!vp->is_vacant())
        {
            cell.push_back(std::make_pair(vp, coord));
        }
        else
        {
            throw NotFound("1");
        }
    }

    void update_matrix(const coordinate_type& from_coord,
        const coordinate_type& to_coord,
        boost::shared_ptr<VoxelPool> vp)
    {
        const matrix_type::size_type from_idx(coordinate2index(from_coord)),
            to_idx(coordinate2index(to_coord));
        if (from_idx == to_idx)
        {
            cell_type& cell(matrix_[from_idx]);
            cell_type::iterator i(find_from_cell(from_coord, cell));
            if (i == cell.end())
            {
                throw NotFound("2");
            }
            (*i).first = vp;
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
            matrix_[to_idx].push_back(std::make_pair(vp, to_coord));
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

    Integer count_voxels(
        const boost::shared_ptr<VoxelPool>& vp) const
    {
        Integer count(0);
        utils::pair_first_element_unary_predicator<boost::shared_ptr<VoxelPool>, coordinate_type> pred(vp);

        for (matrix_type::const_iterator i(matrix_.begin());
            i != matrix_.end(); ++i)
        {
            count += static_cast<Integer>(
                std::count_if((*i).begin(), (*i).end(), pred));
        }
        return count;
    }

    // /**
    //  * Change the Species at v.coordinate() to v.species.
    //  * The ParticleID must be kept after this update.
    //  */
    // virtual void update_voxel(const Voxel& v)
    // {
    //     const coordinate_type coord(v.coordinate());
    //     // VoxelPool* src_vp(get_voxel_pool(coord));
    //     VoxelPool* src_vp(find_voxel_pool(coord));
    //     VoxelPool* new_vp(get_voxel_pool(v));

    //     if (src_vp->with_voxels() != new_vp->with_voxels())
    //     {
    //         throw NotSupported("ParticleID is needed/lost.");
    //     }

    //     new_vp->add_voxel(src_vp->pop(coord));
    //     update_matrix(coord, new_vp);
    // }

    virtual bool update_voxel(const ParticleID& pid, const Voxel& v);

    virtual std::pair<ParticleID, Voxel> get_voxel_at(const coordinate_type& coord) const
    {
        boost::shared_ptr<const VoxelPool> vp(get_voxel_pool_at_(coord));
        return std::make_pair(
            vp->get_particle_id(coord),
            Voxel(vp->species(), coord, vp->radius(), vp->D(), get_location_serial(vp)));
    }

    virtual bool remove_voxel(const ParticleID& pid)
    {
        std::pair<boost::shared_ptr<VoxelPool>, coordinate_type> target(__get_coordinate(pid));
        if (target.second != -1)
        {
            boost::shared_ptr<VoxelPool> vp(target.first);
            const coordinate_type coord(target.second);
            if (!vp->remove_voxel_if_exists(coord))
            {
                return false;
            }

            vp->location()->add_voxel(coordinate_id_pair_type(ParticleID(), coord));
            update_matrix(coord, vp->location());
            return true;
        }
        return false;
    }

    virtual bool remove_voxel(const coordinate_type& coord)
    {
        boost::shared_ptr<VoxelPool> vp(get_voxel_pool_at_(coord));
        if (vp->is_vacant())
        {
            return false;
        }

        if (vp->remove_voxel_if_exists(coord))
        {
            // ???
            update_matrix(coord, vacant_);
            return true;
        }
        return true;
    }

    virtual bool move(
        const coordinate_type& src, const coordinate_type& dest,
        const std::size_t candidate=0)
    {
        coordinate_type tmp_dest(dest);
        if (src == tmp_dest)
        {
            return false;
        }

        boost::shared_ptr<VoxelPool> src_vp(get_voxel_pool_at_(src));
        if (src_vp->is_vacant())
        {
            return true;
        }

        boost::shared_ptr<VoxelPool> dest_vp(get_voxel_pool_at_(tmp_dest));
        if (dest_vp == border_)
        {
            return false;
        }
        else if (dest_vp == periodic_)
        {
            tmp_dest = periodic_transpose(tmp_dest);
            dest_vp = get_voxel_pool_at_(tmp_dest);
        }

        if (dest_vp != src_vp->location())
        {
            return false;
        }

        src_vp->replace_voxel(src, tmp_dest);
        dest_vp->replace_voxel(tmp_dest, src);
        if (!dest_vp->is_vacant())
        {
            update_matrix(src, dest_vp);
            update_matrix(tmp_dest, src_vp);
        }
        else
        {
            update_matrix(src, tmp_dest, src_vp);
        }
        return true;
    }

    virtual bool can_move(const coordinate_type& src, const coordinate_type& dest) const
    {
        if (src == dest)
        {
            return false;
        }

        boost::shared_ptr<const VoxelPool> src_vp(get_voxel_pool_at_(src));
        if (src_vp->is_vacant())
        {
            return false;
        }

        boost::shared_ptr<VoxelPool> dest_vp(get_voxel_pool_at_(dest));

        if (dest_vp == border_)
        {
            return false;
        }

        if (dest_vp == periodic_)
        {
            dest_vp = get_voxel_pool_at_(periodic_transpose(dest));
        }

        return (dest_vp == src_vp->location());
    }

    virtual const Particle particle_at(const coordinate_type& coord) const
    {
        boost::shared_ptr<const VoxelPool> vp(get_voxel_pool_at_(coord));
        return Particle(
            vp->species(), coordinate2position(coord), vp->radius(), vp->D());
    }

    boost::shared_ptr<VoxelPool> get_voxel_pool_at_(const coordinate_type& coord) const;
    VoxelPool* get_voxel_pool_at(const coordinate_type& coord) const
    {
        return get_voxel_pool_at_(coord).get();
    }

    coordinate_type get_neighbor_boundary(
        const coordinate_type& coord, const Integer& nrand) const
    {
        coordinate_type const dest = get_neighbor(coord, nrand);
        return (!is_periodic_ || is_inside(dest) ? dest : periodic_transpose(dest));
    }

    virtual void add_structure(const Species& sp,
        const boost::shared_ptr<const Shape>& s, const std::string loc);
    virtual const boost::shared_ptr<const Shape>& get_structure(const Species& sp) const;
    virtual const Shape::dimension_kind get_structure_dimension(const Species& sp) const;

    virtual std::pair<coordinate_type, bool> move_to_neighbor(
        VoxelPool* const& from_vp, VoxelPool* const& loc,
        coordinate_id_pair_type& info, const Integer nrand);

    virtual Integer num_molecules(const Species& sp) const;

    /**
     */

#ifdef WITH_HDF5
    /*
     * HDF5 Save
     */
    void save_hdf5(H5::Group* root) const
    {
        // save_lattice_space(*this, root, "LatticeSpaceCellListImpl");
        throw NotSupported("LatticeSpaceCellListImpl::save_hdf5 is not supported yet.");
    }

    void load_hdf5(const H5::Group& root)
    {
        // load_lattice_space(root, this);
        // load_lattice_space(root, this, "LatticeSpaceCellListImpl");
        throw NotSupported("LatticeSpaceCellListImpl::load_hdf5 is not supported yet.");
    }
#endif

protected:

    boost::shared_ptr<VoxelPool> get_voxel_pool_(const Voxel& v);

    VoxelPool* get_voxel_pool(const Voxel& v)
    {
        return get_voxel_pool_(v).get();
    }

    std::pair<boost::shared_ptr<VoxelPool>, coordinate_type>
        __get_coordinate(const ParticleID& pid);
    std::pair<boost::shared_ptr<const VoxelPool>, coordinate_type>
        __get_coordinate(const ParticleID& pid) const;

    bool make_molecular_type(
        const Species& sp, Real radius, Real D, const std::string loc);
    bool make_structure_type(
        const Species& sp, Shape::dimension_kind dimension, const std::string loc);

protected:

    bool is_periodic_;

    boost::shared_ptr<VoxelPool> border_;
    boost::shared_ptr<VoxelPool> periodic_;

    Integer3 matrix_sizes_, cell_sizes_;
    matrix_type matrix_;
    structure_container_type structures_;
};

} // ecell4

#endif /* ECELL4_LATTICE_SPACE_CELL_LIST_IMPL_HPP */
