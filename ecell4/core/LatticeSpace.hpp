#ifndef __ECELL4_LATTICE_SPACE_HPP
#define __ECELL4_LATTICE_SPACE_HPP

#include <vector>
#include <set>
#include <map>
#include <stdexcept>

#include "Space.hpp"
#include "Integer3.hpp"
#include "LatticeSpaceHDF5Writer.hpp"
#include "MolecularTypeBase.hpp"
#include "MolecularType.hpp"


namespace ecell4
{

class LatticeSpace
    : public Space
{
public:

    typedef MolecularTypeBase::particle_info particle_info_type;
    typedef MolecularTypeBase::private_coordinate_type private_coordinate_type;
    typedef private_coordinate_type coordinate_type;

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

    const Real& t() const
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

    /**
      */

    Real voxel_radius() const
    {
        return voxel_radius_;
    }

    Real voxel_volume() const
    {
        const Real r(voxel_radius_);
        return 4.0 * sqrt(2.0) * r * r * r;
    }

    virtual const Integer col_size() const = 0;
    virtual const Integer row_size() const = 0;
    virtual const Integer layer_size() const = 0;

    virtual std::vector<Species> list_species() const = 0;

    virtual Integer num_voxels_exact(const Species& sp) const = 0;
    virtual Integer num_voxels(const Species& sp) const = 0;
    virtual Integer num_voxels() const = 0;
    virtual bool has_voxel(const ParticleID& pid) const = 0;

    virtual std::vector<std::pair<ParticleID, Voxel> > list_voxels() const = 0;
    virtual std::vector<std::pair<ParticleID, Voxel> >
        list_voxels(const Species& sp) const = 0;
    virtual std::vector<std::pair<ParticleID, Voxel> >
        list_voxels_exact(const Species& sp) const = 0;

    virtual bool update_voxel_private(const Voxel& v) = 0;
    virtual bool update_voxel_private(const ParticleID& pid, const Voxel& v) = 0;
    virtual std::pair<ParticleID, Voxel> get_voxel(const ParticleID& pid) const = 0;
    virtual bool remove_voxel(const ParticleID& pid) = 0;
    virtual bool remove_voxel_private(const private_coordinate_type& coord) = 0;
    virtual bool move(const coordinate_type& from, const coordinate_type& to) = 0;
    virtual const Particle particle_at(const coordinate_type& coord) const = 0;

    virtual MolecularTypeBase* find_molecular_type(const Species& sp) = 0;
    virtual MolecularTypeBase* get_molecular_type(
        const private_coordinate_type& coord) = 0;

    virtual private_coordinate_type get_neighbor(
        const private_coordinate_type& private_coord, const Integer& nrand) const = 0;

    /**
     Coordinate transformations
     */

    virtual coordinate_type global2coord(const Integer3& global) const = 0;
    virtual coordinate_type global2private_coord(const Integer3& global) const = 0;
    virtual Integer3 private_coord2global(const private_coordinate_type& coord) const = 0;
    virtual Real3 private2position(
        const private_coordinate_type& private_coord) const = 0;
    virtual private_coordinate_type position2private(const Real3& pos) const = 0;
    virtual Real3 coordinate2position(const coordinate_type& coord) const = 0;
    virtual coordinate_type position2coordinate(const Real3& pos) const = 0;
    virtual private_coordinate_type coord2private(const coordinate_type& cood) const = 0;
    virtual coordinate_type private2coord(
        const private_coordinate_type& private_coord) const = 0;
    virtual Integer3 coord2global(coordinate_type coord) const = 0;
    virtual Real3 global2position(const Integer3& global) const = 0;
    virtual Integer3 position2global(const Real3& pos) const = 0;

    virtual bool on_structure(const Voxel& v) = 0;

    virtual std::pair<private_coordinate_type, bool> move_to_neighbor(
        MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
        particle_info_type& info, const Integer nrand) = 0;

    /**
      */

    virtual Integer num_molecules() const
    {
        return num_voxels();
    }

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

    virtual bool update_voxel(const ParticleID& pid, const Voxel& v)
    {
        return update_voxel_private(pid,
            Voxel(v.species(), coord2private(v.coordinate()),
                v.radius(), v.D(), v.loc()));
    }

    virtual bool update_particle(const ParticleID& pid, const Particle& p)
    {
        return update_voxel_private(pid, Voxel(p.species(),
            position2private(p.position()), p.radius(), p.D()));
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

    inline Integer size() const
    {
        return col_size() * row_size() * layer_size();
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
        const Real3& edge_lengths, const Real& voxel_radius)
        : base_type(voxel_radius), edge_lengths_(edge_lengths)
    {
        //XXX: derived from SpatiocyteStepper::setLatticeProperties()
        HCP_L = voxel_radius_ / sqrt(3.0);
        HCP_X = voxel_radius_ * sqrt(8.0 / 3.0); // Lx
        HCP_Y = voxel_radius_ * sqrt(3.0); // Ly

        const Real lengthX = edge_lengths_[0];
        const Real lengthY = edge_lengths_[1];
        const Real lengthZ = edge_lengths_[2];

        row_size_ = (Integer)rint((lengthZ / 2) / voxel_radius_) + 2;
        layer_size_ = (Integer)rint(lengthY / HCP_Y) + 2;
        col_size_ = (Integer)rint(lengthX / HCP_X) + 2;
    }

    virtual ~LatticeSpaceBase()
    {
        ; // do nothing
    }

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

    /**
     Coordinate transformations
     */

    static coordinate_type __global2coord(
        const Integer3& global,
        const Integer& num_col, const Integer& num_row, const Integer& num_layer)
    {
        return global.row + num_row * (global.col + num_col * global.layer);
    }

    static Integer3 __coord2global(
        const coordinate_type& coord,
        const Integer& num_col, const Integer& num_row, const Integer& num_layer)
    {
        const Integer NUM_COLROW(num_row * num_col);
        const Integer LAYER(coord / NUM_COLROW);
        const Integer SURPLUS(coord - LAYER * NUM_COLROW);
        const Integer COL(SURPLUS / num_row);
        const Integer3 retval(COL, SURPLUS - COL * num_row, LAYER);
        return retval;
    }

    /** global -> */

    coordinate_type global2coord(const Integer3& global) const
    {
        return __global2coord(global, col_size(), row_size(), layer_size());
    }

    coordinate_type global2private_coord(const Integer3& global) const
    {
        const Integer3 g(global.col + 1, global.row + 1, global.layer + 1);
        return __global2coord(global, col_size_, row_size_, layer_size_);
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

    /** -> global */

    Integer3 coord2global(coordinate_type coord) const
    {
        return __coord2global(coord, col_size(), row_size(), layer_size());
    }

    Integer3 private_coord2global(const private_coordinate_type& coord) const
    {
        const Integer3 private_global(
            __coord2global(coord, col_size_, row_size_, layer_size_));
        const Integer3 retval(
            private_global.col - 1, private_global.row - 1, private_global.layer - 1);
        return retval;
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

    /** others */

    Real3 coordinate2position(const coordinate_type& coord) const
    {
        // return global2position(private_coord2global(
        //     global2private_coord(coord2global(coord))));
        return global2position(coord2global(coord));
    }

    private_coordinate_type coord2private(const coordinate_type& coord) const
    {
        return global2private_coord(coord2global(coord));
    }

    Real3 private2position(
        const private_coordinate_type& private_coord) const
    {
        return global2position(private_coord2global(private_coord));
    }

    private_coordinate_type position2private(const Real3& pos) const
    {
        return global2private_coord(position2global(pos));
    }

    coordinate_type position2coordinate(const Real3& pos) const
    {
        return global2coord(position2global(pos));
    }

    coordinate_type private2coord(
        const private_coordinate_type& private_coord) const
    {
        return global2coord(private_coord2global(private_coord));
    }

    private_coordinate_type get_neighbor(
        const private_coordinate_type& private_coord, const Integer& nrand) const
    {
        const Integer NUM_COLROW(col_size_ * row_size_);
        const Integer NUM_ROW(row_size_);
        const bool odd_col(((private_coord % NUM_COLROW) / NUM_ROW) & 1);
        const bool odd_lay((private_coord / NUM_COLROW) & 1);

        switch (nrand)
        {
        case 1:
            return private_coord + 1;
        case 2:
            return private_coord + (odd_col ^ odd_lay) - NUM_ROW - 1;
        case 3:
            return private_coord + (odd_col ^ odd_lay) - NUM_ROW;
        case 4:
            return private_coord + (odd_col ^ odd_lay) + NUM_ROW - 1;
        case 5:
            return private_coord + (odd_col ^ odd_lay) + NUM_ROW;
        case 6:
            return private_coord - (2 * odd_col - 1) * NUM_COLROW - NUM_ROW;
        case 7:
            return private_coord - (2 * odd_col - 1) * NUM_COLROW + NUM_ROW;
        case 8:
            return private_coord + (odd_col ^ odd_lay) - NUM_COLROW - 1;
        case 9:
            return private_coord + (odd_col ^ odd_lay) - NUM_COLROW;
        case 10:
            return private_coord + (odd_col ^ odd_lay) + NUM_COLROW - 1;
        case 11:
            return private_coord + (odd_col ^ odd_lay) + NUM_COLROW;
        }
        return private_coord - 1; // nrand == 0
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
    typedef base_type::private_coordinate_type private_coordinate_type;
    typedef base_type::private_coordinate_type coordinate_type;

    typedef std::map<Species, MolecularType> spmap;
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
    virtual bool remove_voxel_private(const private_coordinate_type& coord);

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

    virtual Integer num_voxels_exact(const Species& sp) const;
    virtual Integer num_voxels(const Species& sp) const;
    virtual Integer num_voxels() const;

    virtual Integer num_molecules() const
    {
        return LatticeSpace::num_molecules();
    }

    virtual Integer num_molecules(const Species& sp) const; //XXX:

    // bool update_voxel(const ParticleID& pid, const Voxel& v);
    virtual bool update_voxel_private(const Voxel& v);
    virtual bool update_voxel_private(const ParticleID& pid, const Voxel& v);

    std::vector<Species> list_species() const;
    const Species& find_species(std::string name) const;
    std::vector<coordinate_type> list_coords(const Species& sp) const;
    std::vector<coordinate_type> list_coords_exact(const Species& sp) const;
    virtual MolecularTypeBase* find_molecular_type(const Species& sp);
    // MolecularTypeBase* find_molecular_type(const std::string name);
    virtual MolecularTypeBase* get_molecular_type(const private_coordinate_type& coord);
    // bool register_species(const Species& sp);
    // bool update_molecule(private_coordinate_type coord, const Species& species);
    // bool add_molecule(const Species& sp, private_coordinate_type coord, const ParticleID& pid);
    virtual bool move(const coordinate_type& from, const coordinate_type& to);

    std::pair<private_coordinate_type, bool> move_to_neighbor(
        private_coordinate_type coord, Integer nrand);
    std::pair<private_coordinate_type, bool> move_to_neighbor(
        particle_info_type& info, Integer nrand);
    std::pair<private_coordinate_type, bool> move_to_neighbor(
        MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
        particle_info_type& info, const Integer nrand);

    inline bool is_periodic() const
    {
        return is_periodic_;
    }

    /*
     * HDF5 Save
     */
    void save(H5::Group* root) const
    {
        save_lattice_space(*this, root);
    }

    void load(const H5::Group& root)
    {
        load_lattice_space(root, this);
    }

    void reset(const Real3& edge_lengths, const Real& voxel_radius,
        const bool is_periodic)
    {
        edge_lengths_ = edge_lengths;
        voxel_radius_ = voxel_radius;
        is_periodic_ = is_periodic;

        voxels_.clear();
        spmap_.clear();
        set_lattice_properties(is_periodic_);
    }

    /*
     * Coordinate transformations
     */
    // virtual Integer3 private_coord2global(const private_coordinate_type& coord) const;
    // virtual coordinate_type global2private_coord(const Integer3& global) const;
    // virtual private_coordinate_type coord2private(const coordinate_type& cood) const;
    // virtual coordinate_type private2coord(
    //     const private_coordinate_type& private_coord) const;
    // virtual Real3 coordinate2position(const coordinate_type& coord) const;
    // virtual coordinate_type position2coordinate(const Real3& pos) const;
    // virtual Real3 private2position(const private_coordinate_type& private_coord) const;
    // virtual private_coordinate_type position2private(const Real3& pos) const;
    // virtual coordinate_type global2coord(const Integer3& global) const;
    // virtual Integer3 coord2global(coordinate_type coord) const;
    // virtual Real3 global2position(const Integer3& global) const;
    // virtual Integer3 position2global(const Real3& pos) const;

    // std::vector<private_coordinate_type> get_neighbors(
    //         private_coordinate_type coord) const;
    // virtual private_coordinate_type get_neighbor(
    //     private_coordinate_type private_coord, Integer nrand) const;

    virtual const Particle particle_at(const coordinate_type& coord) const
    {
        return particle_at_private(coord2private(coord));
    }

public:

    /*
     * Coordinate transformations
     */
    private_coordinate_type global2coord_(const Integer3& global,
            Integer col_size, Integer row_size, Integer layer_size) const;
    const Integer3 coord2global_(coordinate_type coord,
            Integer col_size, Integer row_size, Integer layer_size) const;

    const Integer3 private_coord2private_global(
            const private_coordinate_type privatre_coord) const;

    private_coordinate_type apply_boundary_(
            const private_coordinate_type& private_coord) const;


    const spmap& molecular_types() const
    {
        return spmap_;
    }

    bool is_inside(private_coordinate_type coord) const;

    virtual bool on_structure(const Voxel& v);

protected:

    std::pair<spmap::iterator, bool> __get_molecular_type(const Voxel& v);
    MolecularTypeBase* get_molecular_type(const Voxel& v);

    void set_lattice_properties(const bool is_periodic);
    std::pair<private_coordinate_type, bool> move_(
            private_coordinate_type private_from, private_coordinate_type private_to);
    std::pair<private_coordinate_type, bool> move_(
            particle_info_type& info, private_coordinate_type private_to);
    private_coordinate_type get_coord(const ParticleID& pid) const;
    const Particle particle_at_private(private_coordinate_type coord) const;
    bool is_in_range(coordinate_type coord) const;
    bool is_in_range_private(private_coordinate_type coord) const;

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
