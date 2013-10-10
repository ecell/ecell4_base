#ifndef __ECELL4_LATTICE_SPACE_HPP
#define __ECELL4_LATTICE_SPACE_HPP

#include "Space.hpp"
#include "MolecularType.hpp"
#include "SParticle.hpp"
#include "SerialIDGenerator.hpp"
#include <vector>
#include <set>
#include <map>

namespace ecell4
{

//Lattice type:
#define HCP_LATTICE   0
#define CUBIC_LATTICE 1

//Comp dimensions:
#define VOLUME  3
#define SURFACE 2
#define LINE    1

//Comp geometries:
#define CUBOID        0
#define ELLIPSOID     1
#define CYLINDER      2
#define ROD           3
#define PYRAMID       4
#define ERYTHROCYTE   5

//CUBOID Comp surface boundary conditions:
#define REFLECTIVE     0
#define PERIODIC       1
#define UNIPERIODIC    2
#define REMOVE_UPPER   3
#define REMOVE_LOWER   4
#define REMOVE_BOTH    5

//The 12 adjoining voxels of a voxel in the HCP lattice:
#define NORTH    0
#define SOUTH    1
#define NW       2
#define SW       3
#define NE       4
#define SE       5
#define EAST     6
#define WEST     7
#define DORSALN  8
#define DORSALS  9
#define VENTRALN 10
#define VENTRALS 11

//The 6 adjoining voxels of a voxel in the CUBIC lattice:
#define DORSAL   4
#define VENTRAL  5

#define INNER     0
#define OUTER     1
#define IMMED     2
#define EXTEND    3
#define SHARED    4

// Polymerization parameters
#define LARGE_DISTANCE 50
#define MAX_MONOMER_OVERLAP 0.2
#define MAX_IMMEDIATE_DISTANCE 0.2
#define BIG_NUMBER 1e+20

class LatticeSpace
    : public Space
{
protected:

    typedef std::map<ParticleID, Voxel> lattice_container_type;
    typedef std::vector<MolecularType> molecular_type_set;

public:

    LatticeSpace();
    ~LatticeSpace();
    /*
     * methods of Space class
     */
    Integer num_species() const;
    bool has_species(const Species& sp) const;
    Integer num_molecules(const Species& sp)const;
    const Position3& edge_lengths() const;
    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;
    bool has_particle(const ParticleID& pid) const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles(const Species& sp) const;


    /*
     * original methods
     */
    bool update_sparticle(const ParticleID pid, const SParticle spcl);
    void remove_sparticle(ParticleID pid);
    Species add_molecular_type(const std::string name);
    MolecularType& get_molecular_type(Species& sp);
    void coord2global(Integer coord, Integer& global_row,
            Integer& global_layer, Integer& global_col) const;
    const Position3 coord2position(Integer coord) const;

protected:

    Voxel& voxel_as(ParticleID pid);
    Voxel& voxel_at(Integer coord);
    Integer position2coord(Integer row, Integer layer, Integer col);
    const Particle voxel2particle(const Voxel& voxel) const
    {
        const MolecularType* ptr_mt = voxel.ptr_mt;
        const Species& sp = ptr_mt->species();
        const Position3& pos = coord2position(voxel.coord);
        const Real& radius = 0;
        const Real& D = 0;
        Particle particle(sp, pos, radius, D);
        return particle;
    }
    void exchange_coords(Voxel& voxel0, Voxel& voxel1);

    /*
     * Spatiocyte methods
     */
    void set_lattice_properties();
    void concatenate_voxel(Voxel& voxel,
            Integer row, Integer layer, Integer col);
    void concatenate_rows(Voxel& voxel, Integer coord,
            Integer row, Integer layer, Integer col);
    void concatenate_layers(Voxel& voxel, Integer coord,
            Integer row, Integer layer, Integer col);
    void concatenate_cols(Voxel& voxel, Integer coord,
            Integer row, Integer layer, Integer col);

protected:

    Real theNormalizedVoxelRadius;
    Real theHCPl, theHCPx, theHCPy;
    Position3 theCenterPoint;

    Integer lattice_type_;
    lattice_container_type lattice_;
    molecular_type_set molecular_types_;
    Integer theNullCoord;

    MolecularType VACANT_TYPE_;

    Integer adjoining_size_;    // theAdjoiningCoordSize
    Position3 edge_lengths_;
    Integer row_size_, layer_size_, col_size_;
    Integer stride_;

};

}

#endif
