#ifndef __ECELL4_LATTICE_SPACE_HPP
#define __ECELL4_LATTICE_SPACE_HPP

#include <vector>
#include <set>
#include <map>
#include <stdexcept>

#include "Space.hpp"
#include "Integer3.hpp"
#ifdef HDF5
#include "LatticeSpaceHDF5Writer.hpp"
#endif
#include "MolecularTypeBase.hpp"
#include "MolecularType.hpp"


namespace ecell4
{

class LatticeSpace
    : public Space
{
public:

    typedef MolecularTypeBase::particle_info particle_info;
    typedef MolecularTypeBase::private_coordinate_type private_coordinate_type;
    typedef private_coordinate_type coordinate_type;

    typedef std::map<Species, MolecularType> spmap;

    // typedef Integer coordinate_type;
    // typedef coordinate_type private_coordinate_type;
    // typedef std::pair<private_coordinate_type, ParticleID> particle_info;


    typedef std::vector<MolecularTypeBase*> voxel_container;

public:

    LatticeSpace(
        const Real3& edge_lengths, const Real& voxel_radius,
        const bool is_periodic = true);
    ~LatticeSpace();

    /*
     * Space APIs
     *
     * using ParticleID, Species and Posision3
     */

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

    const Real3& edge_lengths() const;

    const Real volume() const
    {
        return edge_lengths_[0] * edge_lengths_[1] * edge_lengths_[2];
    }

    Integer num_species() const;
    Integer num_molecules()const;
    Integer num_molecules(const Species& sp)const;
    Integer num_molecules_exact(const Species& sp)const;
    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;
    Integer num_particles_exact(const Species& sp) const;

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
    bool has_particle(const ParticleID& pid) const;
    bool has_voxel(const ParticleID& pid) const;
    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const;

    std::vector<std::pair<ParticleID, Particle> >
        list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles_exact(const Species& sp) const;

    bool remove_particle(const ParticleID& pid);
    bool remove_voxel(const ParticleID& pid);
    bool remove_voxel_private(const private_coordinate_type coord);

    bool update_particle(const ParticleID& pid, const Particle& p);
    bool update_structure(const Particle& p);

    /*
     * for Simulator
     *
     * using Species and coordinate_type
     */
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels(const Species& sp) const;
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels_exact(const Species& sp) const;
    std::pair<ParticleID, Voxel> get_voxel(const ParticleID& pid) const;

    Integer num_voxels_exact(const Species& sp) const;
    Integer num_voxels(const Species& sp) const;
    Integer num_voxels() const;

    bool update_voxel(const ParticleID& pid, const Voxel& v);
    bool update_voxel_private(const Voxel& v);
    bool update_voxel_private(const ParticleID& pid, const Voxel& v);

    std::vector<Species> list_species() const;
    const Species& find_species(std::string name) const;
    std::vector<coordinate_type> list_coords(const Species& sp) const;
    std::vector<coordinate_type> list_coords_exact(const Species& sp) const;
    MolecularTypeBase* find_molecular_type(const Species& sp);
    // MolecularTypeBase* find_molecular_type(const std::string name);
    MolecularTypeBase* get_molecular_type(private_coordinate_type coord) const;
    // bool register_species(const Species& sp);
    // bool update_molecule(private_coordinate_type coord, const Species& species);
    // bool add_molecule(const Species& sp, private_coordinate_type coord, const ParticleID& pid);
    bool move(coordinate_type from, coordinate_type to);

    std::pair<private_coordinate_type, bool> move_to_neighbor(
        private_coordinate_type coord, Integer nrand);
    std::pair<private_coordinate_type, bool> move_to_neighbor(
        particle_info& info, Integer nrand);
    std::pair<private_coordinate_type, bool> move_to_neighbor(
        MolecularTypeBase* const& from_mt, MolecularTypeBase* const& loc,
        particle_info& info, const Integer nrand);

    inline bool is_periodic() const
    {
        return is_periodic_;
    }

    Real voxel_radius() const
    {
        return voxel_radius_;
    }

    inline Integer col_size() const
    {
        return col_size_ - 2;
    }

    inline Integer row_size() const
    {
        return row_size_ - 2;
    }

    inline Integer layer_size() const
    {
        return layer_size_ - 2;
    }

    inline Integer size() const
    {
        return col_size() * row_size() * layer_size();
    }

#ifdef HDF5
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
#endif

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
    coordinate_type global2coord(const Integer3& global) const;
    coordinate_type global2private_coord(const Integer3& global) const;
    const Integer3 coord2global(coordinate_type coord) const;
    const Integer3 private_coord2global(private_coordinate_type coord) const;

    private_coordinate_type coord2private(coordinate_type cood) const;
    coordinate_type private2coord(private_coordinate_type private_coord) const;

    const Real3 coordinate2position(coordinate_type coord) const;
    coordinate_type position2coordinate(const Real3& pos) const;

    const Real3 private2position(private_coordinate_type private_coord) const;
    private_coordinate_type position2private(const Real3& pos) const;

    std::vector<private_coordinate_type> get_neighbors(
            private_coordinate_type coord) const;
    private_coordinate_type get_neighbor(
            private_coordinate_type private_coord, Integer nrand) const;

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

    const Real3 global2position(const Integer3& global) const;
    const Integer3 position2global(const Real3& pos) const;

    private_coordinate_type apply_boundary_(
            const private_coordinate_type& private_coord) const;


    const spmap& molecular_types() const
    {
        return spmap_;
    }

    bool is_inside(private_coordinate_type coord) const;

    bool on_structure(const Voxel& v);

protected:

    std::pair<spmap::iterator, bool> __get_molecular_type(const Voxel& v);
    MolecularTypeBase* get_molecular_type(const Voxel& v);

    void set_lattice_properties(const bool is_periodic);
    std::pair<private_coordinate_type, bool> move_(
            private_coordinate_type private_from, private_coordinate_type private_to);
    std::pair<private_coordinate_type, bool> move_(
            particle_info& info, private_coordinate_type private_to);
    private_coordinate_type get_coord(const ParticleID& pid) const;
    const Particle particle_at(private_coordinate_type coord) const;
    bool is_in_range(coordinate_type coord) const;
    bool is_in_range_private(private_coordinate_type coord) const;

protected:

    Real voxel_radius_;
    Real3 edge_lengths_;
    Real t_;
    bool is_periodic_;

    Real HCP_L, HCP_X, HCP_Y;

    spmap spmap_;
    voxel_container voxels_;

    MolecularTypeBase* vacant_;
    MolecularTypeBase* border_;
    MolecularTypeBase* periodic_;

    Integer row_size_, layer_size_, col_size_;
};

}

#endif
