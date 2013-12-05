#ifndef __ECELL4_LATTICE_SPACE_HPP
#define __ECELL4_LATTICE_SPACE_HPP

#include "Space.hpp"
#include "MolecularType.hpp"
#include "VacantType.hpp"
#include "SParticle.hpp"
#include "Global.hpp"
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

    typedef std::map<Species, MolecularType> spmap;
    typedef std::vector<MolecularType*> voxel_container;


public:

    LatticeSpace();
    LatticeSpace(const Position3& edge_lengths);
    ~LatticeSpace();

    /*
     * APIs
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
     * using Species, SParticle and Coord
     */
    std::vector<SParticle> list_sparticles() const;
    std::vector<SParticle> list_sparticles(const Species& sp) const;

    bool move(Integer from, Integer to);
    bool has_particle(Integer coord) const;


protected:

    void set_lattice_properties();

    Integer get_coord(const ParticleID& pid) const;
    MolecularType* voxel_at(Integer coord) const;
    MolecularType* get_molecular_type(const Species& sp);

    /*
     * coordinate transformations
     */
    const Global coord2global(Integer coord) const;
    const Position3 coord2position(Integer coord) const;

    Integer global2coord(const Global& global) const;
    const Position3 global2position(const Global& global) const;

    Integer position2coord(const Position3& pos) const;
    const Global position2global(const Position3& pos) const;

    const Particle particle_at(Integer coord) const;


protected:

    Real theNormalizedVoxelRadius;
    Real theHCPl, theHCPx, theHCPy;
    Position3 theCenterPoint;

    Integer lattice_type_;
    spmap spmap_;
    voxel_container voxels_;

    Position3 edge_lengths_;
    Integer row_size_, layer_size_, col_size_;
    Integer stride_;

};

}

#endif
