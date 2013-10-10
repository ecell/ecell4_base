#include "LatticeSpace.hpp"

namespace ecell4
{

LatticeSpace::LatticeSpace()
    : VACANT_TYPE_("VACANT")
{
    set_lattice_properties();
}

LatticeSpace::~LatticeSpace()
{
}

Integer LatticeSpace::num_species() const
{
    return molecular_types_.size();
}

bool LatticeSpace::has_species(const Species& sp) const
{
    for (molecular_type_set::const_iterator itr(molecular_types_.begin());
            itr != molecular_types_.end(); ++itr)
    {
        const MolecularType& mt(*itr);
        if (mt.species() == sp)
        {
            return true;
        }
    }
    return false;
}

Integer LatticeSpace::num_molecules(const Species& sp) const
{
    return num_particles(sp) * 1;
}

const Position3& LatticeSpace::edge_lengths() const
{
    return edge_lengths_;
}

Integer LatticeSpace::num_particles() const
{
    Integer count(0);
    for (lattice_container_type::const_iterator i(lattice_.begin());
            i != lattice_.end(); ++i)
    {
        const Voxel voxel((*i).second);
        const MolecularType* p_mt(voxel.ptr_mt);
        if (p_mt != &VACANT_TYPE_)
        {
            count++;
        }
    }
    return count;
}

Integer LatticeSpace::num_particles(const Species& sp) const
{
    for (molecular_type_set::const_iterator i(molecular_types_.begin());
            i != molecular_types_.end(); ++i)
    {
        const MolecularType& mt = *i;
        if (mt.species() == sp)
            return mt.voxels().size();
    }
    return 0;
}

bool LatticeSpace::has_particle(const ParticleID& pid) const
{
    return lattice_.find(pid) != lattice_.end();
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles() const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    for (lattice_container_type::const_iterator i(lattice_.begin());
            i != lattice_.end(); ++i)
    {
        const Voxel voxel((*i).second);
        retval.push_back(std::pair<ParticleID, Particle>(
                    voxel.id, voxel2particle(voxel)));
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    for (molecular_type_set::const_iterator i(molecular_types_.begin());
            i != molecular_types_.end(); ++i)
    {
        const MolecularType& mt = *i;
        if (mt.species() == sp)
        {
            const MolecularType::voxel_container_type vct = mt.voxels();
            for (MolecularType::voxel_container_type::const_iterator itr(vct.begin());
                    itr != vct.end(); ++itr)

            {
                const Voxel& voxel = (*itr).second;
                retval.push_back(std::pair<ParticleID, Particle>(
                            voxel.id, voxel2particle(voxel)));
            }
        }
    }
    return retval;
}

/*
 * Spatiocyte methods
 */

/*
 * derived from SpatiocyteStepper::setLatticeProperties()
 */
void LatticeSpace::set_lattice_properties()
{
    adjoining_size_ = 12;
    theHCPl = theNormalizedVoxelRadius/sqrt(3);
    theHCPx = theNormalizedVoxelRadius*sqrt(8.0/3); //Lx
    theHCPy = theNormalizedVoxelRadius*sqrt(3); //Ly

    Integer lengthX = edge_lengths_[0];
    Integer lengthY = edge_lengths_[1];
    Integer lengthZ = edge_lengths_[2];

    theCenterPoint[2] = lengthZ / 2 + 4 *
        theNormalizedVoxelRadius; //row
    theCenterPoint[1] = lengthY / 2 + 2 * theHCPy; //layer
    theCenterPoint[0] = lengthX / 2 + 2 * theHCPx; //column

    row_size_ = (Integer)rint((theCenterPoint[2])/
                              (theNormalizedVoxelRadius));
    layer_size_ = (Integer)rint((theCenterPoint[1]*2)/theHCPy);
    col_size_ = (Integer)rint((theCenterPoint[0]*2)/theHCPx);

    Integer coord(0);
    SerialIDGenerator<ParticleID> sidgen;
    for (Integer coord(0); coord < row_size_ * layer_size_ * col_size_; ++coord)
    {
        ParticleID pid(sidgen());
        Voxel voxel(pid, coord, NULL);
        VACANT_TYPE_.addVoxel(voxel);
        lattice_.insert(lattice_container_type::value_type(pid, voxel));
    }
    for (Integer coord(0); coord < row_size_ * layer_size_ * col_size_; ++coord)
    {
        Voxel& voxel(voxel_at(coord));
        Position3 pos(coord2position(coord));
        concatenate_voxel(voxel, pos[0], pos[1], pos[2]);
    }
    theNullCoord = row_size_ * layer_size_ * col_size_;
}

void LatticeSpace::concatenate_voxel(Voxel& voxel,
        Integer row, Integer layer, Integer col)
{
    Integer coord(position2coord(row, layer, col));
    if (row > 0)
    {
        concatenate_rows(voxel, coord, row-1, layer, col);
    }
    if (layer > 0)
    {
        concatenate_layers(voxel, coord, row, layer-1, col);
    }
    if (col > 0)
    {
        concatenate_cols(voxel, coord, row, layer, col-1);
    }
}

void LatticeSpace::concatenate_rows(Voxel& voxel, Integer coord,
        Integer row, Integer layer, Integer col)
{
    Integer b(position2coord(row, layer, col));
    Voxel& adjoining_voxel(voxel_at(b));
    voxel.adjoiningVoxels[NORTH] = &adjoining_voxel;
    adjoining_voxel.adjoiningVoxels[SOUTH] = &voxel;
}

void LatticeSpace::concatenate_layers(Voxel& voxel, Integer coord,
        Integer row, Integer layer, Integer col)
{
    Integer b(position2coord(row, layer, col));
    Voxel& adjoining_voxel(voxel_at(b));
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            if ((layer + 1) % 2 + col % 2 == 1)
            {
                voxel.setAdjoiningVoxel(VENTRALN, &adjoining_voxel);
                adjoining_voxel.setAdjoiningVoxel(DORSALS, &voxel);
                if (row > row_size_ - 1)
                {
                    Integer c(position2coord(row+1, layer, col));
                    Voxel& adjoining_voxel(voxel_at(c));
                    voxel.setAdjoiningVoxel(VENTRALS, &adjoining_voxel);
                    adjoining_voxel.setAdjoiningVoxel(DORSALN, &voxel);
                }
            }
            else
            {
                voxel.setAdjoiningVoxel(VENTRALS, &adjoining_voxel);
                adjoining_voxel.setAdjoiningVoxel(DORSALN, &voxel);
                if (row > 0)
                {
                    Integer c(position2coord(row-1, layer, col));
                    Voxel& adjoining_voxel(voxel_at(c));
                    voxel.setAdjoiningVoxel(VENTRALN, &adjoining_voxel);
                    adjoining_voxel.setAdjoiningVoxel(DORSALS, &voxel);
                }
            }
            break;
        case CUBIC_LATTICE:
            voxel.setAdjoiningVoxel(VENTRAL, &adjoining_voxel);
            adjoining_voxel.setAdjoiningVoxel(DORSAL, &voxel);
            break;
    }
}

void LatticeSpace::concatenate_cols(Voxel& voxel, Integer coord,
        Integer row, Integer layer, Integer col)
{
    Integer b(position2coord(row, layer, col));
    Voxel& adjoining_voxel(voxel_at(b));
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            if (layer % 2 == 0)
            {
                if ((col + 1) % 2 == 1)
                {
                    voxel.setAdjoiningVoxel(NW, &adjoining_voxel);
                    adjoining_voxel.setAdjoiningVoxel(SE, &voxel);
                    if (row < row_size_ - 1)
                    {
                        Integer c(position2coord(row+1, layer, col));
                        Voxel& adjoining_voxel(voxel_at(c));
                        voxel.setAdjoiningVoxel(SW, &adjoining_voxel);
                        adjoining_voxel.setAdjoiningVoxel(NE, &voxel);
                    }
                    if (layer < layer_size_ - 1)
                    {
                        Integer c(position2coord(row, layer+1, col));
                        Voxel& adjoining_voxel(voxel_at(c));
                        voxel.setAdjoiningVoxel(WEST, &adjoining_voxel);
                        adjoining_voxel.setAdjoiningVoxel(EAST, &voxel);
                    }
                }
                else
                {
                    voxel.adjoiningVoxels[SW] = &adjoining_voxel;
                    adjoining_voxel.adjoiningVoxels[NE] = &voxel;
                    if (row > 0)
                    {
                        Integer c(position2coord(row-1, layer, col));
                        Voxel& adjoining_voxel(voxel_at(c));
                        voxel.setAdjoiningVoxel(NW, &adjoining_voxel);
                        adjoining_voxel.setAdjoiningVoxel(SE, &voxel);
                    }
                    if (layer > 0)
                    {
                        Integer c(position2coord(row, layer-1, col));
                        Voxel& adjoining_voxel(voxel_at(c));
                        voxel.setAdjoiningVoxel(WEST, &adjoining_voxel);
                        adjoining_voxel.setAdjoiningVoxel(EAST, &voxel);
                    }
                }
            }
            else
            {
                if ((col + 1) % 2 == 1)
                {
                    voxel.setAdjoiningVoxel(SW, &adjoining_voxel);
                    adjoining_voxel.setAdjoiningVoxel(NE, &voxel);
                    if (row > 0)
                    {
                        Integer c(position2coord(row-1, layer, col));
                        Voxel& adjoining_voxel(voxel_at(c));
                        voxel.setAdjoiningVoxel(NW, &adjoining_voxel);
                        adjoining_voxel.setAdjoiningVoxel(SE, &voxel);
                    }
                    if (layer < layer_size_ - 1)
                    {
                        Integer c(position2coord(row, layer+1, col));
                        Voxel& adjoining_voxel(voxel_at(c));
                        voxel.setAdjoiningVoxel(WEST, &adjoining_voxel);
                        adjoining_voxel.setAdjoiningVoxel(EAST, &voxel);
                    }
                }
                else
                {
                    voxel.setAdjoiningVoxel(NW, &adjoining_voxel);
                    adjoining_voxel.setAdjoiningVoxel(SE, &voxel);
                    if (row < row_size_ -1)
                    {
                        Integer c(position2coord(row+1, layer, col));
                        Voxel& adjoining_voxel(voxel_at(c));
                        voxel.setAdjoiningVoxel(SW, &adjoining_voxel);
                        adjoining_voxel.setAdjoiningVoxel(NE, &voxel);
                    }
                    if (layer > 0)
                    {
                        Integer c(position2coord(row, layer-1, col));
                        Voxel& adjoining_voxel(voxel_at(c));
                        voxel.setAdjoiningVoxel(WEST, &adjoining_voxel);
                        adjoining_voxel.setAdjoiningVoxel(EAST, &voxel);
                    }
                }
            }
            break;
        case CUBIC_LATTICE:
            voxel.setAdjoiningVoxel(WEST, &adjoining_voxel);
            adjoining_voxel.setAdjoiningVoxel(EAST, &voxel);
    }
}

/*
 * original methods
 */

bool LatticeSpace::update_sparticle(ParticleID pid, SParticle spcl)
{
    Voxel& src = voxel_as(pid);
    Voxel& dest = voxel_at(spcl.coord);
    if (dest.ptr_mt != &VACANT_TYPE_)
    {
        return false;
    }

    MolecularType* src_ptr_mt(src.ptr_mt);
    MolecularType& dest_mt = get_molecular_type(spcl.species);

    src_ptr_mt->removeVoxel(src.id);
    dest_mt.addVoxel(src);

    exchange_coords(src, dest);

    return true;
}

void LatticeSpace::exchange_coords(Voxel& voxel0, Voxel& voxel1)
{
    Integer tmp_coord = voxel0.coord;
    Integer tmp_diffuse_size = voxel0.diffuse_size;
    Voxel** tmp_adjoiningVoxels = voxel0.adjoiningVoxels;

    voxel0.coord = voxel1.coord;
    voxel0.diffuse_size = voxel1.diffuse_size;
    voxel0.adjoiningVoxels = voxel1.adjoiningVoxels;

    voxel1.coord = tmp_coord;
    voxel1.diffuse_size = tmp_diffuse_size;
    voxel1.adjoiningVoxels = tmp_adjoiningVoxels;

    Position3 pos0 = coord2position(voxel0.coord);
    Position3 pos1 = coord2position(voxel1.coord);
    concatenate_voxel(voxel0, pos0[0], pos0[1], pos0[2]);
    concatenate_voxel(voxel1, pos1[0], pos1[1], pos1[2]);
}

void LatticeSpace::remove_sparticle(ParticleID pid)
{
    Voxel& voxel = voxel_as(pid);
    MolecularType* p_mt(voxel.ptr_mt);
    p_mt->removeVoxel(voxel.id);
    VACANT_TYPE_.addVoxel(voxel);
}

Species LatticeSpace::add_molecular_type(const std::string name)
{
    const MolecularType mt(name);
    molecular_types_.push_back(mt);
    return mt.species();
}

MolecularType& LatticeSpace::get_molecular_type(Species& sp)
{
    for (molecular_type_set::iterator i(molecular_types_.begin());
            i != molecular_types_.end(); ++i)
    {
        const MolecularType& mt = (*i);
        if (mt.species() == sp)
        {
            return const_cast<MolecularType&>(mt);
        }
    }
    throw "Exception in get_molerular_type";
}

void LatticeSpace::coord2global(Integer aCoord, Integer& aGlobalRow,
        Integer& aGlobalLayer, Integer& aGlobalCol) const
{
    aGlobalCol = aCoord / (row_size_ * layer_size_);
    aGlobalLayer = (aCoord % (row_size_ * layer_size_)) / row_size_;
    aGlobalRow = (aCoord % (row_size_ * layer_size_)) % row_size_;
}

const Position3 LatticeSpace::coord2position(Integer coord) const
{
    Integer aGlobalCol;
    Integer aGlobalLayer;
    Integer aGlobalRow;
    coord2global(coord, aGlobalRow, aGlobalLayer, aGlobalCol);
    //the center point of a voxel
    Position3 position;
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            position[1] = (aGlobalCol % 2) * theHCPl + theHCPy * aGlobalLayer;
            position[2] = aGlobalRow * 2 * theNormalizedVoxelRadius +
            ((aGlobalLayer + aGlobalCol) % 2) * theNormalizedVoxelRadius;
            position[0] = aGlobalCol * theHCPx;
            break;
        case CUBIC_LATTICE:
            position[1] = aGlobalLayer * 2 * theNormalizedVoxelRadius;
            position[2] = aGlobalRow * 2 * theNormalizedVoxelRadius;
            position[0] = aGlobalCol * 2 * theNormalizedVoxelRadius;
            break;
    }
    return position;
}

/*
* Protected methods
*/
Voxel& LatticeSpace::voxel_as(ParticleID pid)
{
    lattice_container_type::iterator i(lattice_.find(pid));
    if (i == lattice_.end())
    {
        throw "Exception: Not in lattice_";
    }
    return (*i).second;
}

Voxel& LatticeSpace::voxel_at(Integer coord)
{
    for (lattice_container_type::iterator i(this->lattice_.begin());
            i != lattice_.end(); ++i)
    {
        if ((*i).second.coord = coord)
        {
            return (*i).second;
        }
    }
}

Integer LatticeSpace::position2coord(Integer row, Integer layer, Integer col)
{
    return row +
        row_size_ * layer +
        row_size_ * layer_size_ * col;
}

}
