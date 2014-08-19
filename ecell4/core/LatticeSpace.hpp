#ifndef __ECELL4_LATTICE_SPACE_HPP
#define __ECELL4_LATTICE_SPACE_HPP

#include <vector>
#include <set>
#include <map>
#include <stdexcept>

#include "Space.hpp"
#include "Global.hpp"
#include "LatticeSpaceHDF5Writer.hpp"


namespace ecell4
{

class MolecularTypeBase;
class MolecularType;

class LatticeSpace
    : public Space
{
protected:

    typedef std::vector<MolecularTypeBase*> voxel_container;

public:

    typedef std::map<Species, MolecularType> spmap;
    typedef Integer coordinate_type;
    typedef coordinate_type private_coordinate_type;

    typedef std::pair<private_coordinate_type, ParticleID> particle_info;

    LatticeSpace(
        const Position3& edge_lengths, const Real& voxel_radius,
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

    const Position3& edge_lengths() const;

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
    std::vector<coordinate_type> list_coords(const Species& sp) const;
    std::vector<coordinate_type> list_coords_exact(const Species& sp) const;
    MolecularTypeBase* find_molecular_type(const Species& sp);
    MolecularTypeBase* get_molecular_type(private_coordinate_type coord) const;
    // bool register_species(const Species& sp);
    // bool update_molecule(private_coordinate_type coord, const Species& species);
    // bool add_molecule(const Species& sp, private_coordinate_type coord, const ParticleID& pid);
    bool move(coordinate_type from, coordinate_type to);
    std::pair<private_coordinate_type, bool> move_to_neighbor(private_coordinate_type coord, Integer nrand);
    std::pair<private_coordinate_type, bool> move_to_neighbor(particle_info& info, Integer nrand);

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

    void cleanup(const Position3& edge_lengths, const Real& voxel_radius,
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
    coordinate_type global2coord(const Global& global) const;
    coordinate_type global2private_coord(const Global& global) const;
    const Global coord2global(coordinate_type coord) const;
    const Global private_coord2global(private_coordinate_type coord) const;

    private_coordinate_type coord2private(coordinate_type cood) const;
    coordinate_type private2coord(private_coordinate_type private_coord) const;

    const Position3 coordinate2position(coordinate_type coord) const;
    coordinate_type position2coordinate(const Position3& pos) const;

    private_coordinate_type get_neighbor(
            private_coordinate_type private_coord, Integer nrand) const;

public:

    /*
     * Coordinate transformations
     */
    private_coordinate_type global2coord_(const Global& global,
            Integer col_size, Integer row_size, Integer layer_size) const;
    const Global coord2global_(coordinate_type coord,
            Integer col_size, Integer row_size, Integer layer_size) const;

    const Global private_coord2private_global(
            const private_coordinate_type privatre_coord) const;
    const private_coordinate_type private_global2private_coord(
            const Global private_global) const;

    const Position3 global2position(const Global& global) const;
    const Global position2global(const Position3& pos) const;

    private_coordinate_type apply_boundary_(
            const private_coordinate_type& private_coord) const;

    private_coordinate_type position2private_coord(const Position3& pos) const;

    const spmap& molecular_types() const
    {
        return spmap_;
    }

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
    Position3 edge_lengths_;
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
