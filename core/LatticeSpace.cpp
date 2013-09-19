#include "LatticeSpace.hpp"

namespace ecell4
{

LatticeSpace::LatticeSpace()
    : VACANT_TYPE_("VACANT")
{
    set_lattice_properties();
}

Integer LatticeSpace::num_species() const
{
    const species_pointer_set sps = get_species_set();
    return sps.size();
}

bool LatticeSpace::has_species(const Species& sp) const
{
    const species_pointer_set sps = get_species_set();
    return (sps.find(const_cast<Species*>(&sp)) != sps.end());
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
    Integer count = 0;
    for (Integer i = 0; i < lattice_.size(); ++i)
    {
        const Voxel& voxel = lattice_.at(i);
        const MolecularType* p_mt = voxel.p_molecule_type;
        if (p_mt != &VACANT_TYPE_)
        {
            count++;
        }
    }
    return count;
}

Integer LatticeSpace::num_particles(const Species& sp) const
{
    for (molecular_type_set::iterator i(molecular_types_.begin());
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
    for (Integer i = 0; i < lattice_.size(); ++i)
    {
        const Voxel& voxel = lattice_.at(i);
        if (voxel.id == pid)
            return true;
    }
    return false;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles() const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    for (Integer i = 0; i < lattice_.size(); ++i)
    {
        const Voxel& voxel = lattice_.at(i);
        retval.push_back(std::pair<ParticleID, Particle>(
                    voxel.id, voxel2particle(voxel)));
    }
    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    for (molecular_type_set::iterator i(molecular_types_.begin());
            i != molecular_types_.end(); ++i)
    {
        const MolecularType& mt = *i;
        if (mt.species() == sp)
        {
            const MolecularType::voxel_container_type vct = mt.voxels();
            for (Integer i = 0; i < vct.size(); ++i)
            {
                const Voxel& voxel = *(vct.at(i));
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

    lattice_.resize(row_size_ * layer_size_ * col_size_ + 1);
    for (lattice_container_type::iterator i(lattice_.begin());
            i != lattice_.end(); ++i)
    {
        (*i).p_molecule_type = &VACANT_TYPE_;
    }
    theNullCoord = row_size_ * layer_size_ * col_size_;
}

/*
 * original methods
 */

bool LatticeSpace::update_sparticle(ParticleID pid, SParticle spcl)
{
    Voxel& src = voxel_as(pid);
    Voxel& dest = voxel_at(spcl.coord);
    if (dest.p_molecule_type != &VACANT_TYPE_)
    {
        return false;
    }
    MolecularType* src_mt_p = src.p_molecule_type;

    src_mt_p->removeVoxel(src);
    MolecularType& mt = get_molecular_type(spcl.species);
    mt.addVoxel(&src);

    dest.p_molecule_type = src_mt_p;
    dest.id = src.id;

    src.p_molecule_type = &VACANT_TYPE_;
    src.id = VACANT_ID_;
    return true;
}

void LatticeSpace::remove_sparticle(ParticleID pid)
{
    Voxel& v = voxel_as(pid);
    MolecularType* p_mt = v.p_molecule_type;
    p_mt->removeVoxel(v);
    v.p_molecule_type = &VACANT_TYPE_;
    v.id = VACANT_ID_;
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

const LatticeSpace::species_pointer_set LatticeSpace::get_species_set() const
{
    species_pointer_set retval;
    for (molecular_type_set::iterator i(molecular_types_.begin());
            i != molecular_types_.end(); ++i)
    {
        const MolecularType& mt = *i;
        if (&mt != &VACANT_TYPE_)
        {
            Species s = (*i).species();
            retval.insert(&s);
        }
    }
    return retval;
}

Voxel& LatticeSpace::voxel_as(ParticleID pid)
{
    for (lattice_container_type::iterator i(this->lattice_.begin());
            i != lattice_.end(); ++i)
    {
        if ((*i).id == pid)
        {
            return *i;
        }
    }
}

Voxel& LatticeSpace::voxel_at(Integer coord)
{
    for (lattice_container_type::iterator i(this->lattice_.begin());
            i != lattice_.end(); ++i)
    {
        if ((*i).coord = coord)
        {
            return *i;
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
