#ifndef __ECELL4_LATTICE_SPACE_HPP
#define __ECELL4_LATTICE_SPACE_HPP

#include <vector>
#include <set>
#include <map>
#include <stdexcept>
#include "Space.hpp"
#include "Global.hpp"
#include "Voxel.hpp"
#include "LatticeSpaceHDF5Writer.hpp"

namespace ecell4
{

class MolecularTypeBase;
class MolecularType;

class LatticeSpace
    : public Space
{
protected:

    typedef std::map<Species, MolecularType> spmap;
    typedef std::vector<MolecularTypeBase*> voxel_container;

public:

    typedef Integer coordinate_type;
    typedef coordinate_type private_coordinate_type;

    typedef std::pair<private_coordinate_type, ParticleID> particle_info;

    LatticeSpace(const Position3& edge_lengths, const Real& voxel_radius, const bool is_periodic);
    ~LatticeSpace();

    /*
     * Space APIs
     *
     * using ParticleID, Species and Posision3
     */
    const Position3& edge_lengths() const;

    Integer num_species() const;
    Integer num_molecules(const Species& sp)const;
    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;

    bool has_species(const Species& sp) const;
    bool has_particle(const ParticleID& pid) const;

    std::vector<std::pair<ParticleID, Particle> >
        list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles(const Species& sp) const;

    bool update_particle(const ParticleID& pid, const Particle& p);

    /*
     * for Simulator
     *
     * using Species and coordinate_type
     */
    std::vector<std::pair<ParticleID, Voxel> >
        list_voxels(const Species& sp) const;

    Integer num_voxels(const Species& sp) const
    {
        spmap::const_iterator itr(spmap_.find(sp));
        if (itr == spmap_.end())
        {
            return 0;
        }
        const MolecularTypeBase* mt(&((*itr).second));
        return mt->size();
    }

    Integer num_voxels() const
    {
        Integer retval(0);
        for (spmap::const_iterator itr(spmap_.begin()); itr != spmap_.end(); ++itr)
        {
            retval += (*itr).second.size();
        }
        return retval;
    }

    bool update_voxel(const ParticleID& pid, const Voxel& v)
    {
        const Coord to_coord(inner2general(v.coordinate()));
        if (to_coord < 0 && to_coord >= row_size_ * col_size_ * layer_size_) //XXX: is_in_range
        {
            return false;
        }

        MolecularTypeBase* dest_mt(get_molecular_type(v.species()));
        const Coord from_coord(get_coord(pid));
        if (from_coord != -1)
        {
            MolecularTypeBase* src_mt(voxels_.at(from_coord));
            src_mt->removeVoxel(from_coord);
            voxel_container::iterator itr(voxels_.begin() + from_coord);
            voxels_.erase(itr);
            voxels_.insert(itr, vacant_);
        }

        dest_mt->addVoxel(MolecularTypeBase::particle_info(to_coord, pid));
        voxel_container::iterator itr(voxels_.begin() + to_coord);
        (*itr) = dest_mt;
        return true;
    }

    std::vector<Species> list_species() const;
    std::vector<coordinate_type> list_coords(const Species& sp) const;
    MolecularTypeBase* get_molecular_type(const Species& sp);
    MolecularTypeBase* get_molecular_type(private_coordinate_type coord) const;
    bool add_species(const Species& sp);
    bool add_molecule(const Species& sp, private_coordinate_type coord, const ParticleID& pid);
    bool remove_molecule(const private_coordinate_type coord);
    bool move(coordinate_type from, coordinate_type to);
    std::pair<private_coordinate_type, bool> move_to_neighbor(private_coordinate_type coord, Integer nrand);
    std::pair<private_coordinate_type, bool> move_to_neighbor(particle_info& info, Integer nrand);
    bool update_molecule(private_coordinate_type coord, const Species& species);

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
     * Coordinate transformations
     */
    coordinate_type global2coord(const Global& global) const;
    coordinate_type global2private_coord(const Global& global) const;
    const Global coord2global(coordinate_type coord) const;
    const Global private_coord2global(private_coordinate_type coord) const;

    private_coordinate_type coord2private(coordinate_type cood) const;
    coordinate_type private2coord(private_coordinate_type private_coord) const;

    /*
     * HDF5 Save
     */
    void save(H5::Group* root) const
    {
        // LatticeSpaceHDF5Writer<LatticeSpace> writer(*this);
        // writer.save(fout, hdf5path);
        save_lattice_space(*this, root);
    }

    void load(const H5::Group& root)
    {
        load_lattice_space(root, this);
    }

    void cleanup(const Position3& edge_lengths, const Real& voxel_radius)
    {
        edge_lengths_ = edge_lengths;
        theNormalizedVoxelRadius = voxel_radius;

        voxels_.clear();
        spmap_.clear();
        set_lattice_properties();
    }

protected:

    void set_lattice_properties(const bool is_periodic);
    private_coordinate_type get_neighbor(
            private_coordinate_type private_coord, Integer nrand) const;
    std::pair<private_coordinate_type, bool> move_(
            private_coordinate_type private_from, private_coordinate_type private_to);
    std::pair<private_coordinate_type, bool> move_(
            particle_info& info, private_coordinate_type private_to);
    private_coordinate_type get_coord(const ParticleID& pid) const;
    const Particle particle_at(private_coordinate_type coord) const;
    bool is_in_range(coordinate_type coord) const;

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

    const Position3 coord2position(coordinate_type coord) const;
    coordinate_type position2coord(const Position3& pos) const;

    private_coordinate_type apply_boundary_(
            const private_coordinate_type& private_coord) const;

protected:

    Real voxel_radius_;
    Position3 edge_lengths_;
    Real t_;

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
