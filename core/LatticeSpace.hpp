#ifndef __ECELL4_LATTICE_SPACE_HPP
#define __ECELL4_LATTICE_SPACE_HPP

#include "Space.hpp"
#include "MolecularType.hpp"
#include "SParticle.hpp"
#include <set>

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

// Simply Comp (temporary)
struct Comp
{
    bool isIntersectParent;
    bool isIntersectRoot;
// unsigned
    Integer dimension;
    Integer vacantID; //remove this
    Integer interfaceID;
// int
    Integer enclosed;
    Integer geometry;
    Integer xyPlane;
    Integer xzPlane;
    Integer yzPlane;
// unsigned
    Integer minRow;
    Integer minCol;
    Integer minLayer;
    Integer maxRow;
    Integer maxCol;
    Integer maxLayer;
// double
    Real lengthX;
    Real lengthY;
    Real lengthZ;
    Real originX;
    Real originY;
    Real originZ;
    Real rotateX;
    Real rotateY;
    Real rotateZ;
    Real specVolume;
    Real specArea;
    Real actualVolume;
    Real actualArea;
//    System* system;
    Comp* surfaceSub;
    //Even if there are many adjacent diffusive compartents, use only one single
    //common id. So there is only one common diffusive Comp:
    Comp* diffusiveComp;
    Position3 centerPoint;
    Species* vacantSpecies;
    std::vector<Comp*> allSubs;
    std::vector<Comp*> immediateSubs;
    std::vector<Comp*> intersectPeers;
    std::vector<Comp*> intersectLowerPeers;
    std::vector<Comp*> lineSubs;
    std::vector<Species*> species;
// std::vector<unsigned int>
    std::vector<Integer> adjoinCount;
};

class LatticeSpace
    : public Space
{
protected:

    typedef std::vector<Voxel> lattice_container_type;
    typedef std::set<MolecularType> molecular_type_set;
    typedef std::set<Species*> species_pointer_set;

public:

    LatticeSpace();
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
     * Spatiocyte methods
     */
    void set_lattice_properties();
    void construct_lattice();
    void concatenate_voxel(Voxel& voxel,
            Integer row, Integer layer, Integer col);
    void concatenate_rows(Voxel& voxel,
            Integer row, Integer layer, Integer col);
    void concatenate_layers(Voxel& voxel,
            Integer row, Integer layer, Integer col);
    void concatenate_cols(Voxel& voxel,
            Integer row, Integer layer, Integer col);
    bool is_inside_coord(Integer coord, Comp* comp, Real delta);

    /*
     * original methods
     */
    bool update_sparticle(ParticleID pid, SParticle spcl);
    void remove_sparticle(ParticleID pid);
    void coord2global(Integer coord, Integer& global_row,
            Integer& global_layer, Integer& global_col) const;
    const Position3 coord2position(Integer coord) const;

protected:

    const species_pointer_set get_species_set() const;
    Voxel& voxel_as(ParticleID pid);
    Voxel& voxel_at(Integer coord);
    Integer position2coord(Integer row, Integer layer, Integer col);
    const Particle voxel2particle(const Voxel& voxel) const
    {
        const MolecularType* mt = voxel.p_molecule_type;
        const Species& sp = mt->species();
        const Position3& pos = coord2position(voxel.coord);
        const Real& radius = 0;
        const Real& D = 0;
        Particle particle(sp, pos, radius, D);
        return particle;
    }

protected:

    Real theNormalizedVoxelRadius;
    Real theHCPl, theHCPx, theHCPy;
    Position3 theCenterPoint;

    Integer lattice_type_;
    lattice_container_type lattice_;
    molecular_type_set molecular_types_;
    Integer theNullCoord;

    MolecularType VACANT_TYPE_;
    ParticleID VACANT_ID_;

    Integer adjoining_size_;    // theAdjoiningCoordSize
    Position3 edge_lengths_;
    Integer row_size_, layer_size_, col_size_;
    Integer stride_;

};

}

#endif
