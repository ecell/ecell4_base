#include "LatticeSpace.hpp"

/*
 * Only Supported HCP_Lattice
 */

namespace ecell4
{

LatticeSpace::LatticeSpace()
{

}

Integer LatticeSpace::num_species()
{
    species_set sp_set = species_set()
    return static_cast<Integer>(sp_set.size());
}

bool has_speceis(const Species& sp) const
{
    species_set sp_set = species_set();
    return sp_set.find(sp) != sp_set.end();
}

const Position3& LatticeSpace::edge_lengths() const
{
    return edge_lengths_;
}

Integer LatticeSpace::num_particles() const
{
    return static_cast<Integer>(lattice_.size());
}

Integer LatticeSpace::num_particles(const Species$ sp) const
{
    Integer count(0);
    for (it = lattice_.begin(); it != lattice_.end(); i++)
    {
        const Voxel& voxel = it;
        const MoleculeType& m_type = voxel.molecule_type;
        const Species& species = m_type.species();
        if (species == sp)
            count++;
    }
    return count;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles() const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    voxel_container_type::iterator it;
    for (it = lattice_.begin(); it != lattice_.end(); i++)
    {
        const Voxel& voxel= *it;
        const MoleculeType& m_type = voxel.molecule_type;
        const Species& species = m_type.species();
        const Position3& pos = coord2position(voxel);
        const Particle p(species, pos, &radius, &D);
        retval.push_back(std::pair<ParticleID, Particle>(it.id, p));
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    voxel_container_type::iterator it;
    for (it = lattice_.begin(); it != lattice_.end(); i++)
    {
        const Voxel& voxel= *it;
        const MoleculeType& m_type = voxel.molecule_type;
        const Species& species = m_type.species();
        if (species == sp)
            continue;
        const Position3& pos = coord2position(voxel);
        const Particle p(species, pos, &radius, &D);
        retval.push_back(std::pair<ParticleID, Particle>(it.id, p));
    }
    return retval;
}

void LatticeSpace::update_sparticle(ParticleID id, SParticle sparticle)
{
    Voxel& src = voxel_as(id);
    Voxel& dest = voxel_at(sparticle.coord);
}

Voxel& LatticeSpace::voxel_at(Integer coord) const
{
    for (voxel_container_type::iterator i(lattice_.begin());
            i != lattice_.end(); ++i)
    {
        if ((*i).coord == coord)
            return (*i);
    }
}

Voxel& LatticeSpace::voxel_as(ParticleID id) const
{
    for (voxel_container_type::iterator i(lattice_.begin());
            i != lattice_.end(); ++i)
    {
        if ((*i).id == id)
            return (*i);
    }
}

void LatticeSpace:::set_properties()
{
    center_point_[2] = edge_lengths_[2] / 2 + 4 * NVR;
    center_point_[1] = edge_lengths_[1] / 2 + 2 * HCPy;
    center_point_[0] = edge_lengths_[0] / 2 + 2 * HCPx;

    row_size_ = (Integer)rint(center_point[2] / NVR);
    layer_size_ = (Integer)rint(center_point[1] * 2 / HCPy);
    col_size_ = (Integer)rint(center_point[0] * 2 / HCPx);

    lattice_.resize(row_size_ * layer_size_ * col_size_ + 1);
    Integer null_coord = row_size_ * layer_size_ * col_size_;
    Integer null_id = 0;
    lattice[null_coord].id = null_id;
}

void LatticeSpace::set_adjoining()
{
    Integer coord(0);
    for (voxel_container_type::iterator i(lattice_.begin());
            i != lattice_.end(); ++i, ++)
    {
        (*i).coord = coord;
        (*i).adjoiningVoxels.resize(adjoining_size);
        Integer col(coord / (low_size_ * layer_size_);
        Integer row((coord % (row_size_ * layer_size_)) / row_size_);
        Integer layer((coord % (row_size_ * layer_size_)) % row_size_);
        concatnate_voxel();
    }
}

void LatticeSpace::concatenate_voxel(Voxel& voxel,
    Integer row, Integer layer, Integer col)
{
    if (row > 0)
    {
        concatenate_rows(voxel, row-1, layer, col);
    }
    if (layer > 0)
    {
        concatenate_rows(voxel, row, layer-1, col);
    }
    if (col > 0)
    {
        concatenate_rows(voxel, row, layer, col-1);
    }
}

void concatenate_rows(Voxel& voxel,
        Integer row, Integer layer, Integer col)
{
    Integer b(position2coord(row, layer, col));
    Voxel& adjioning(lattice_.at(b));
    voxel.adjoiningVoxels[NORTH] = adjoining;
    adjioining.adjoiningVoxels[SOURTH] = voxel;
}

void concatenate_layers(Voxel& voxel,
        Integer row, Integer layer, Integer col)
{
    Integer b(position2coord(row, layer, col));
    Voxel& adjoining(lattice_.at(b));

    Integer index[4];
    Integer sign, border;
    if ((layer + 1) % 2 + col % 2 == 1)
    {
        index[0] = VENTRALN;
        index[1] = DORSALS;
        index[2] = VENTRALS;
        index[3] = DORSALN;
        sign = 1;
        border = row_size_ - 1;
    }
    else
    {
        index[0] = VENTRALS;
        index[1] = DORSALN;
        index[2] = VENTRALN;
        index[3] = DORSALS;
        sign = -1;
        border = 0;
    }

    voxel.adjoiningVoxels[index[0]] = adjoining;
    adjoining.adjoiningVoxels[index[1]] = voxel;
    if (row * sign < border * sign)
    {
        Integer c(position2coord(row + sign, layer, col));
        Voxel& adjoining(lattice_.at(c));
        voxel.adjoiningVoxels[index[2]] = adjoining;
        adjoining.adjoiningVoxels[index[3]] = voxel;
    }
}


void concatenate_cols(Voxel& voxel,
        Integer row, Integer layer, Integer col)
{
    Integer b(position2coord(row, layer, col));
    Voxel& adjoining(lattice_.at(b));

    if(layer%2 == 0)
    {
        if((col+1)%2 == 1)
        {
            voxel.adjoiningVoxels[NW] = voxel;
            adjoining.adjoiningVoxels[SE] = adjoining;

            if(row < row_size_ - 1)
            {
                Integer c(position2coord(row + 1, layer, col);
                Voxel& adjoining(lattice_.at(c));
                voxel.adjoiningVoxels[SW] = adjoining;
                adjoining.adjoiningVoxels[NE] = voxel;
            }
            if(layer < layer_size_-1)
            {
                Integer c(position2coord(row, layer + 1, col);
                Voxel& adjoining(lattice_.at(c));
                voxel.adjoiningVoxels[WEST] = adjoining;
                adjoining.adjoiningVoxels[EAST] = voxel;
            }
        }
        else
        {
            voxel.adjoiningVoxels[SW] = adjoining;
            adjoining.adjoiningVoxels[NE] = voxel;

            if(row > 0)
            {
                Integer c(position2coord(row - 1, layer, col));
                Voxel& adjoining(lattice_.at(c));
                voxel.adjoiningVoxels[NW] = adjoining;
                adjoining.adjoiningVoxels[SE] = voxel;
            }
            if(layer > 0)
            {
                Integer c(position2coord(row, layer - 1, col));
                Voxel& adjoining(lattice_.at(c));
                voxel.adjoiningVoxels[WEST] = adjoining;
                adjoining.adjoiningVoxels[EAST] = voxel;
            }
        }
    }
    else
    {
        if((col+1)%2 == 1)
        {
            voxel.adjoiningVoxels[SW] = adjoining;
            adjoining.adjoiningVoxels[NE] = voxel;

            if(row > 0)
            {
                Integer c(position2coord(row - 1, layer, col));
                Voxel& adjoining(lattice_.at(c));
                voxel.adjoiningVoxels[NW] = adjoining;
                adjoining.adjoiningVoxels[SE] = voxel;
            }
            if(layer < layer_size_-1)
            {
                Integer c(position2coord(row, layer + 1, col));
                Voxel& adjoining(lattice_.at(c));
                voxel.adjoiningVoxels[WEST] = adjoining;
                adjoining.adjoiningVoxels[EAST] = voxel;
            }
        }
        else
        {
            voxel.adjoiningVoxels[NW] = adjoining;
            adjoining.adjoiningVoxels[SE] = voxel;
            if(row < row_size_ - 1)
            {
                Integer c(position2coord(row + 1, layer, col));
                Voxel& adjoining(lattice_[c]);
                voxel.adjoiningVoxels[SW] = adjoining;
                adjoining.adjoiningVoxels[NE] = voxel;
            }
            if(layer > 0)
            {
                Integer c(position2coord(row, layer - 1, col));
                Voxel& adjoining(lattice_[c]);
                voxel.adjoiningVoxels[WEST] = adjoining;
                adjoining.adjoiningVoxels[EAST] = voxel;
            }
        }
    }
}

const species_set LatticeSpace::species_set()
{
    species_set sp_set;
    for (molecular_type_set::iterator i(molecular_types_.begin());
            i != molecular_types_.end(); ++I)
    {
        sp_set.insert((*i).species());
    }
    return sp_set;
}

}
