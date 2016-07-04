#ifndef __ECELL4_LATTICE_SPACE_HPP
#define __ECELL4_LATTICE_SPACE_HPP

#include <vector>
#include <set>
#include <map>
#include <stdexcept>

#include "Shape.hpp"
#include "Space.hpp"
#include "Integer3.hpp"

#ifdef WITH_HDF5
#include "LatticeSpaceHDF5Writer.hpp"
#endif

#include "MolecularTypeBase.hpp"
#include "MolecularType.hpp"
#include "Voxel.hpp"

namespace ecell4
{

/**
  * XXX: Just for the temporal use
  */
template <typename Tspace_>
Integer inner2coordinate(const Tspace_& w, const Integer coord)
{
    const Integer num_row(w.row_size());
    const Integer num_col(w.col_size());
    const Integer num_layer(w.layer_size());

    const Integer NUM_COLROW(num_row * num_col);
    const Integer LAYER(coord / NUM_COLROW);
    const Integer SURPLUS(coord - LAYER * NUM_COLROW);
    const Integer COL(SURPLUS / num_row);
    const Integer3 g(COL, SURPLUS - COL * num_row, LAYER);

    return w.global2coordinate(g);
}

class LatticeSpace
    : public Space
{
public:

    typedef Voxel::coordinate_type coordinate_type;
    typedef MolecularTypeBase::coord_id_pair particle_info_type;

public:

    LatticeSpace(const Real& voxel_radius)
        : t_(0.0), voxel_radius_(voxel_radius)
    {
        ;
    }

    virtual ~LatticeSpace()
    {
        ; // do nothing
    }

    // SpaceTraits

    const Real t() const
    {
        return t_;
    }

    void set_t(const Real& t)
    {
        if (t < 0.0)
        {
            throw std::invalid_argument("the time must be positive.");
        }
        t_ = t;
    }

    Real voxel_radius() const
    {
        return voxel_radius_;
    }

    Real voxel_volume() const
    {
        const Real r(voxel_radius_);
        return 4.0 * sqrt(2.0) * r * r * r;
    }

    Real unit_area() const
    {
        const Real r(voxel_radius_);
        return 2.0 * sqrt(3.0) * r * r;
    }

    Real get_volume() const
    {
        return inner_size() * voxel_volume();
    }

    virtual Real3 actual_lengths() const = 0;

    virtual void save(const std::string& filename) const
    {
        throw NotSupported(
            "save(const std::string) is not supported by this space class");
    }

#ifdef WITH_HDF5
    virtual void save_hdf5(H5::Group* root) const
    {
        throw NotSupported(
            "load(H5::Group* root) is not supported by this space class");
    }

    virtual void load_hdf5(const H5::Group& root)
    {
        throw NotSupported(
            "load(const H5::Group& root) is not supported by this space class");
    }
#endif

    /**
      */

    virtual std::vector<Species> list_species() const = 0;

    virtual Integer num_voxels_exact(const Species& sp) const = 0;
    virtual Integer num_voxels(const Species& sp) const = 0;
    virtual Integer num_voxels() const = 0;
    virtual bool has_voxel(const ParticleID& pid) const = 0;

    virtual std::vector<std::pair<ParticleID, Voxel> >
        list_voxels() const = 0;
    virtual std::vector<std::pair<ParticleID, Voxel> >
        list_voxels(const Species& sp) const = 0;
    virtual std::vector<std::pair<ParticleID, Voxel> >
        list_voxels_exact(const Species& sp) const = 0;

    virtual void update_voxel(const Voxel& v) = 0;
    virtual bool update_voxel(const ParticleID& pid, const Voxel& v) = 0;
    virtual bool update_voxel_without_checking(const ParticleID& pid, const Voxel& v)
    {
        throw NotSupported(
            "update_voxel_without_chekcing(const ParticleID&, const Voxel&) is not supported by this space class");
    }

    virtual std::pair<ParticleID, Voxel> get_voxel(const ParticleID& pid) const = 0;
    virtual std::pair<ParticleID, Voxel> get_voxel(const coordinate_type& coord) const = 0;
    virtual bool remove_voxel(const ParticleID& pid) = 0;
    virtual bool remove_voxel(const coordinate_type& coord) = 0;
    virtual bool move(
        const coordinate_type& src, const coordinate_type& dest,
        const std::size_t candidate=0) = 0;
    virtual bool can_move(const coordinate_type& src, const coordinate_type& dest) const;
    virtual const Particle particle_at(const coordinate_type& coord) const = 0;

    virtual MolecularTypeBase* find_molecular_type(const Species& sp) = 0;
    virtual const MolecularTypeBase* find_molecular_type(const Species& sp) const = 0;
    virtual MolecularTypeBase* get_molecular_type(
        const coordinate_type& coord) = 0;
    virtual bool make_structure_type(const Species& sp,
        Shape::dimension_kind dimension, const std::string loc);
    virtual bool make_interface_type(const Species& sp,
        Shape::dimension_kind dimension, const std::string loc);

    virtual bool on_structure(const Voxel& v) = 0;

    virtual std::pair<coordinate_type, bool> move_to_neighbor(
        MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
        particle_info_type& info, const Integer nrand) = 0;

    /**
     Coordinate transformations: See LatticeSpaceBase for the implementation
     */

    virtual const Integer col_size() const = 0;
    virtual const Integer row_size() const = 0;
    virtual const Integer layer_size() const = 0;

    virtual coordinate_type global2coordinate(const Integer3& global) const = 0;
    virtual Integer3 coordinate2global(const coordinate_type& coord) const = 0;
    virtual Real3 coordinate2position(const coordinate_type& coord) const = 0;
    virtual coordinate_type position2coordinate(const Real3& pos) const = 0;
    virtual Real3 global2position(const Integer3& global) const = 0;
    virtual Integer3 position2global(const Real3& pos) const = 0;

    virtual coordinate_type get_neighbor(
        const coordinate_type& coord, const Integer& nrand) const = 0;
    virtual coordinate_type get_neighbor_boundary(
        const coordinate_type& coord, const Integer& nrand) const = 0;

    /**
      */

    virtual Integer num_molecules(const Species& sp) const = 0; //XXX:

    virtual Integer num_molecules_exact(const Species& sp) const
    {
        return num_voxels_exact(sp);
    }

    Integer num_particles() const
    {
        return num_voxels();
    }

    Integer num_particles(const Species& sp) const
    {
        return num_voxels(sp);
    }

    Integer num_particles_exact(const Species& sp) const
    {
        return num_voxels_exact(sp);
    }

    bool has_particle(const ParticleID& pid) const
    {
        return has_voxel(pid);
    }

    virtual std::vector<std::pair<ParticleID, Particle> > list_particles() const
    {
        const std::vector<std::pair<ParticleID, Voxel> > voxels(list_voxels());

        std::vector<std::pair<ParticleID, Particle> > retval;
        retval.reserve(voxels.size());
        for (std::vector<std::pair<ParticleID, Voxel> >::const_iterator
            i(voxels.begin()); i != voxels.end(); ++i)
        {
            const ParticleID& pid((*i).first);
            const Particle p(particle_at((*i).second.coordinate()));
            retval.push_back(std::make_pair(pid, p));
        }
        return retval;
    }

    virtual std::vector<std::pair<ParticleID, Particle> >
        list_particles(const Species& sp) const
    {
        const std::vector<std::pair<ParticleID, Voxel> > voxels(list_voxels(sp));

        std::vector<std::pair<ParticleID, Particle> > retval;
        retval.reserve(voxels.size());
        for (std::vector<std::pair<ParticleID, Voxel> >::const_iterator
            i(voxels.begin()); i != voxels.end(); ++i)
        {
            const ParticleID& pid((*i).first);
            const Particle p(particle_at((*i).second.coordinate()));
            retval.push_back(std::make_pair(pid, p));
        }
        return retval;
    }

    virtual std::vector<std::pair<ParticleID, Particle> >
        list_particles_exact(const Species& sp) const
    {
        const std::vector<std::pair<ParticleID, Voxel> >
            voxels(list_voxels_exact(sp));

        std::vector<std::pair<ParticleID, Particle> > retval;
        retval.reserve(voxels.size());
        for (std::vector<std::pair<ParticleID, Voxel> >::const_iterator
            i(voxels.begin()); i != voxels.end(); ++i)
        {
            const ParticleID& pid((*i).first);
            const Particle p(particle_at((*i).second.coordinate()));
            retval.push_back(std::make_pair(pid, p));
        }
        return retval;
    }

    virtual std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const
    {
        const Voxel v(get_voxel(pid).second);
        return std::make_pair(pid, Particle(
            v.species(), coordinate2position(v.coordinate()), v.radius(), v.D()));
    }

    virtual bool remove_particle(const ParticleID& pid)
    {
        return remove_voxel(pid);
    }

    virtual Integer size() const = 0;
    virtual Integer3 shape() const = 0;

    inline Integer inner_size() const
    {
        return col_size() * row_size() * layer_size();
    }

    inline Integer3 inner_shape() const
    {
        return Integer3(col_size(), row_size(), layer_size());
    }

protected:

    Real t_;
    Real voxel_radius_;
};

class LatticeSpaceBase
    : public LatticeSpace
{
public:

    typedef LatticeSpace base_type;

public:

    LatticeSpaceBase(
        const Real3& edge_lengths, const Real& voxel_radius, const bool is_periodic)
        : base_type(voxel_radius), edge_lengths_(edge_lengths)
    {
        set_lattice_properties(is_periodic);
    }

    virtual ~LatticeSpaceBase()
    {
        ; // do nothing
    }

    virtual void reset(const Real3& edge_lengths, const Real& voxel_radius,
        const bool is_periodic)
    {
        edge_lengths_ = edge_lengths;
        voxel_radius_ = voxel_radius;

        set_lattice_properties(is_periodic);
    }

    void set_lattice_properties(const bool is_periodic);

    /**
     * Primitives
     */

    const Real3& edge_lengths() const
    {
        return edge_lengths_;
    }

    const Real volume() const
    {
        return edge_lengths_[0] * edge_lengths_[1] * edge_lengths_[2];
    }

    virtual const Integer col_size() const
    {
        return col_size_ - 2;
    }

    virtual const Integer row_size() const
    {
        return row_size_ - 2;
    }

    virtual const Integer layer_size() const
    {
        return layer_size_ - 2;
    }

    virtual Real3 actual_lengths() const
    {
        // return Real3(
        //     col_size() * HCP_X,
        //     layer_size() * HCP_Y + (col_size() > 1 ? HCP_L: 0.0),
        //     row_size() * voxel_radius() * 2
        //         + (col_size() > 1 || layer_size() > 1 ? voxel_radius() : 0.0));
        return Real3(
            col_size() * HCP_X, layer_size() * HCP_Y, row_size() * voxel_radius() * 2);
        // const Real sigma(voxel_radius() * 2);
        // return Real3(
        //     (col_size() - 1) * HCP_X + sigma, (layer_size() - 1) * HCP_Y + sigma, row_size() * sigma);
    }

    /**
     Coordinate transformations
     */

    coordinate_type global2coordinate(const Integer3& global) const
    {
        const Integer3 g(global.col + 1, global.row + 1, global.layer + 1);
        return __global2coord(g, col_size_, row_size_, layer_size_);
    }

    Integer3 coordinate2global(const coordinate_type& coord) const
    {
        const Integer3 global(
            __coord2global(coord, col_size_, row_size_, layer_size_));
        const Integer3 retval(
            global.col - 1, global.row - 1, global.layer - 1);
        return retval;
    }

    Real3 coordinate2position(const coordinate_type& coord) const
    {
        return global2position(coordinate2global(coord));
    }

    coordinate_type position2coordinate(const Real3& pos) const
    {
        return global2coordinate(position2global(pos));
    }

    Real3 global2position(const Integer3& global) const
    {
        // the center point of a voxel
        const Real3 pos(
            global.col * HCP_X,
            (global.col % 2) * HCP_L + HCP_Y * global.layer,
            (global.row * 2 + (global.layer + global.col) % 2)
                * voxel_radius_);
        return pos;
    }

    Integer3 position2global(const Real3& pos) const
    {
        const Integer col(round(pos[0] / HCP_X));
        const Integer layer(round((pos[1] - (col % 2) * HCP_L) / HCP_Y));
        const Integer row(round(
            (pos[2] / voxel_radius_ - ((layer + col) % 2)) / 2));
        const Integer3 global(col, row, layer);
        return global;
    }

    coordinate_type get_neighbor(
        const coordinate_type& coord, const Integer& nrand) const
    {
        const Integer NUM_COLROW(col_size_ * row_size_);
        const Integer NUM_ROW(row_size_);
        const bool odd_col(((coord % NUM_COLROW) / NUM_ROW) & 1);
        const bool odd_lay((coord / NUM_COLROW) & 1);

        switch (nrand)
        {
        case 1:
            return coord + 1;
        case 2:
            return coord + (odd_col ^ odd_lay) - NUM_ROW - 1;
        case 3:
            return coord + (odd_col ^ odd_lay) - NUM_ROW;
        case 4:
            return coord + (odd_col ^ odd_lay) + NUM_ROW - 1;
        case 5:
            return coord + (odd_col ^ odd_lay) + NUM_ROW;
        case 6:
            return coord - (2 * odd_col - 1) * NUM_COLROW - NUM_ROW;
        case 7:
            return coord - (2 * odd_col - 1) * NUM_COLROW + NUM_ROW;
        case 8:
            return coord + (odd_col ^ odd_lay) - NUM_COLROW - 1;
        case 9:
            return coord + (odd_col ^ odd_lay) - NUM_COLROW;
        case 10:
            return coord + (odd_col ^ odd_lay) + NUM_COLROW - 1;
        case 11:
            return coord + (odd_col ^ odd_lay) + NUM_COLROW;
        }
        return coord - 1; // nrand == 0
    }

    coordinate_type periodic_transpose(
        const coordinate_type& coord) const
    {
        Integer3 global(coordinate2global(coord));

        global.col = global.col % col_size();
        global.row = global.row % row_size();
        global.layer = global.layer % layer_size();

        global.col = global.col < 0 ? global.col + col_size() : global.col;
        global.row = global.row < 0 ? global.row + row_size() : global.row;
        global.layer = global.layer < 0 ? global.layer + layer_size() : global.layer;

        return global2coordinate(global);
    }

public:

    bool is_in_range(const coordinate_type& coord) const
    {
        return coord >= 0 && coord < row_size_ * col_size_ * layer_size_;
    }

    bool is_inside(const coordinate_type& coord) const
    {
        const Integer3 global(coordinate2global(coord));
        return global.col >= 0 && global.col < col_size()
            && global.row >= 0 && global.row < row_size()
            && global.layer >= 0 && global.layer < layer_size();
    }

    virtual Integer size() const
    {
        return row_size_ * col_size_ * layer_size_;
    }

    virtual Integer3 shape() const
    {
        return Integer3(col_size_, row_size_, layer_size_);
    }

protected:

    static inline Integer __global2coord(
        const Integer3& global,
        const Integer& num_col, const Integer& num_row, const Integer& num_layer)
    {
        return global.row + num_row * (global.col + num_col * global.layer);
    }

    static inline Integer3 __coord2global(
        const Integer& coord,
        const Integer& num_col, const Integer& num_row, const Integer& num_layer)
    {
        const Integer NUM_COLROW(num_row * num_col);
        const Integer LAYER(coord / NUM_COLROW);
        const Integer SURPLUS(coord - LAYER * NUM_COLROW);
        const Integer COL(SURPLUS / num_row);
        const Integer3 retval(COL, SURPLUS - COL * num_row, LAYER);
        return retval;
    }

protected:

    Real3 edge_lengths_;
    Real HCP_L, HCP_X, HCP_Y;
    Integer row_size_, layer_size_, col_size_;
};

class LatticeSpaceVectorImpl
    : public LatticeSpaceBase
{
public:

    typedef LatticeSpaceBase base_type;

    typedef base_type::particle_info_type particle_info_type;
    typedef base_type::coordinate_type coordinate_type;

    typedef std::map<Species, boost::shared_ptr<MolecularType> > spmap;
    typedef std::vector<MolecularTypeBase*> voxel_container;

public:

    LatticeSpaceVectorImpl(
        const Real3& edge_lengths, const Real& voxel_radius,
        const bool is_periodic = true);
    ~LatticeSpaceVectorImpl();

    /*
     * Space APIs
     *
     * using ParticleID, Species and Posision3
     */

    Integer num_species() const;

    virtual Real get_value(const Species& sp) const
    {
        return static_cast<Real>(num_molecules(sp));
    }

    virtual Real get_value_exact(const Species& sp) const
    {
        return static_cast<Real>(num_molecules_exact(sp));
    }

    bool has_species(const Species& sp) const;
    // bool has_species_exact(const Species& sp) const;
    virtual bool has_voxel(const ParticleID& pid) const;

    virtual bool remove_voxel(const ParticleID& pid);
    virtual bool remove_voxel(const coordinate_type& coord);

    bool update_structure(const Particle& p);

    /*
     * for Simulator
     *
     * using Species and coordinate_type
     */
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels() const;
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels(const Species& sp) const;
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels_exact(const Species& sp) const;

    virtual std::pair<ParticleID, Voxel> get_voxel(const ParticleID& pid) const;
    virtual std::pair<ParticleID, Voxel> get_voxel(const coordinate_type& coord) const;

    virtual Integer num_voxels_exact(const Species& sp) const;
    virtual Integer num_voxels(const Species& sp) const;
    virtual Integer num_voxels() const;

    virtual Integer num_molecules(const Species& sp) const; //XXX:

    virtual void update_voxel(const Voxel& v);
    virtual bool update_voxel(const ParticleID& pid, const Voxel& v);
    virtual bool update_voxel_without_checking(const ParticleID& pid, const Voxel& v);

    bool add_voxels(const Species species, std::vector<std::pair<ParticleID, coordinate_type> > voxels);

    std::vector<Species> list_species() const;
    const Species& find_species(std::string name) const;
    std::vector<coordinate_type> list_coords(const Species& sp) const;
    std::vector<coordinate_type> list_coords_exact(const Species& sp) const;
    virtual MolecularTypeBase* find_molecular_type(const Species& sp);
    virtual const MolecularTypeBase* find_molecular_type(const Species& sp) const;
    // MolecularTypeBase* find_molecular_type(const std::string name);
    virtual MolecularTypeBase* get_molecular_type(const coordinate_type& coord);
    // bool update_molecule(coordinate_type coord, const Species& species);
    // bool add_molecule(const Species& sp, coordinate_type coord, const ParticleID& pid);
    virtual bool move(
        const coordinate_type& src, const coordinate_type& dest,
        const std::size_t candidate=0);
    virtual bool can_move(const coordinate_type& src, const coordinate_type& dest) const;

    std::pair<coordinate_type, bool> move_to_neighbor(
        coordinate_type coord, Integer nrand);
    std::pair<coordinate_type, bool> move_to_neighbor(
        particle_info_type& info, Integer nrand);
    std::pair<coordinate_type, bool> move_to_neighbor(
        MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
        particle_info_type& info, const Integer nrand);

    coordinate_type get_neighbor_boundary(
        const coordinate_type& coord, const Integer& nrand) const
    {
        coordinate_type const dest = get_neighbor(coord, nrand);
        MolecularTypeBase* dest_mt(voxels_.at(dest));
        return (dest_mt != periodic_ ? dest : periodic_transpose(dest));
    }

    inline bool is_periodic() const
    {
        return is_periodic_;
    }


#ifdef WITH_HDF5
    /*
     * HDF5 Save
     */
    void save_hdf5(H5::Group* root) const
    {
        save_lattice_space(*this, root);
    }

    void load_hdf5(const H5::Group& root)
    {
        load_lattice_space(root, this);
    }
#endif

    void reset(const Real3& edge_lengths, const Real& voxel_radius,
        const bool is_periodic)
    {
        base_type::reset(edge_lengths, voxel_radius, is_periodic);

        is_periodic_ = is_periodic;
        initialize_voxels(is_periodic_);
    }

    virtual const Particle particle_at(const coordinate_type& coord) const;

    coordinate_type apply_boundary_(
        const coordinate_type& coord) const
    {
        return periodic_transpose(coord);
    }

    virtual bool make_structure_type(const Species& sp,
        Shape::dimension_kind dimension, const std::string loc);
    virtual bool make_interface_type(const Species& sp,
        Shape::dimension_kind dimension, const std::string loc);
    bool make_molecular_type(const Species& sp,
        Real radius, Real D, const std::string loc);

    virtual bool on_structure(const Voxel& v);

protected:

    std::pair<spmap::iterator, bool> __get_molecular_type(const Voxel& v);
    MolecularTypeBase* get_molecular_type(const Voxel& v);

    void initialize_voxels(const bool is_periodic);

    std::pair<coordinate_type, bool> move_(
            coordinate_type from, coordinate_type to,
            const std::size_t candidate=0);
    std::pair<coordinate_type, bool> move_(
            particle_info_type& info, coordinate_type to);
    coordinate_type get_coord(const ParticleID& pid) const;

    Integer count_voxels(const boost::shared_ptr<MolecularType>& mt) const;

protected:

    bool is_periodic_;

    spmap spmap_;
    voxel_container voxels_;

    MolecularTypeBase* vacant_;
    MolecularTypeBase* border_;
    MolecularTypeBase* periodic_;
};

}

#endif
