#ifndef __ECELL4_LATTICE_SPACE_HPP
#define __ECELL4_LATTICE_SPACE_HPP

#include "Space.hpp"
#include "Voxel.hpp"
#include "SParticle.hpp"

#include <set>
#include <math.h>

namespace ecell4
{

class LatticeSpace
    : public Space
{

protected:

    typedef std::vector<Voxel> voxel_container_type;
    typedef std::set<MolecularType> molecular_type_set;
    typedef std::set<Species&> species_set;

public:

    LatticeSpace(Position3 edge_lengths);

    /*
     * override from Space
     */
    Integer num_species() const;
    bool has_species(const Species& sp) const;
    const Position3& edge_lengths() const;
    Integer num_particles() const;
    Integer num_particles(const Species$ sp) const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
        list_particles(const Species& sp) const;

    /*
     * original functions
     */
    //void update_sparticle(ParticleID id, SParticle sparticle);
    //void remove_sparticle(ParticleID id);

    //Voxel& voxel_at(Integer coord) const;
    //Voxel& voxel_as(ParticleID id) const;


protected:

    /*
    void set_properties();
    void set_adjoining();
    void concatenate_voxel(Voxel& voxel,
            Integer row, Integer layer, Integer col);
    void concatenate_rows(Voxel& voxel,
            Integer row, Integer layer, Integer col);
    void concatenate_layers(Voxel& voxel,
            Integer row, Integer layer, Integer col);
    void concatenate_cols(Voxel& voxel,
            Integer row, Integer layer, Integer col);
    */
    const species_set species_set() const;
    const Integer position2coord(Integer row,
            Integer layer, Integer col)
    {
        Integer coord(row +
                row_size_ * layer +
                row_size_ * layer_size_ * col);
        return coord;
    }
    const Position3 coord2position(Voxel& voxel) const
    {
        Integer coord = voxel.id,
                gcol = coord / (row_size_ * layer_size_),
                glayer = (coord % (row_size_ * layer_size_)) / row_size_,
                grow = (coord % (row_size_ * layer_size_)) % row_size_;
        Real y = (gcol % 2) *  HCPl + HCPy * glayer,
             z = grow * 2 * NVR + ((glayer + gcol) % 2) * NVR,
             x = gcol * HCPx;
        return Position3(x, y, z);
    }

protected:

    const static Real NVR = 0.5,    // theNormalizedVoxelRadius
          HCPl = NVR / sqrt(2),
          HCPx = NVR * sqrt(8.0/3),
          HCPy = NVR * sqrt(3);
    const static Integer adjoining_size = 12;
    const static Integer
        NORTH(0),
        SOUTH(1),
        NW(2),
        SW(3),
        NE(4),
        SE(5),
        EAST(6),
        WEST(7),
        DORSALN(8),
        DORSALS(9),
        VENTRALN(10),
        VENTRALS(11);
    Integer col_size_, row_size_, layer_size_;
    Real radius;
    Real D;
    voxel_container_type lattice_;
    molecular_type_set molecular_types_;
    Position3 edge_lengths_;
    Position3 center_point_;

};

}

#endif
