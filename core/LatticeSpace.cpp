#include "LatticeSpace.hpp"

namespace ecell4
{

LatticeSpace::LatticeSpace()
    :  theNormalizedVoxelRadius(0.5)
{
    // example
    edge_lengths_[0] = 10;
    edge_lengths_[1] = 5;
    edge_lengths_[2] = 15;

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
        if (voxel.id != VACANT_ID)
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
    for (lattice_container_type::const_iterator i(lattice_.begin());
            i != lattice_.end(); ++i)
    {
        const Voxel voxel((*i).second);
        if (voxel.id == pid)
            return true;
    }
    return false;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles() const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    for (lattice_container_type::const_iterator i(lattice_.begin());
            i != lattice_.end(); ++i)
    {
        const Voxel voxel((*i).second);
        if (voxel.id == VACANT_ID)
        {
            continue;
        }
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
    lattice_type_ = HCP_LATTICE;

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

    for (Integer coord(0); coord < row_size_ * layer_size_ * col_size_; ++coord)
    {
        ParticleID pid; // default pid
        Voxel voxel(pid, coord, NULL);
        lattice_.insert(lattice_container_type::value_type(coord, voxel));
    }
    //theNullCoord = row_size_ * layer_size_ * col_size_;
}

/*
 * original methods
 */

bool LatticeSpace::update_sparticle(ParticleID pid, SParticle spcl)
{
    Voxel& dest = voxel_at(spcl.coord);
    if (dest.id != VACANT_ID)
    {
        return false;
    }
    MolecularType& dest_mt = get_molecular_type(spcl.species);

    if (has_particle(pid))
    {
        Voxel& src = voxel_as(pid);
        MolecularTypeBase* src_ptr_mt(src.ptr_mt);
        src_ptr_mt->removeVoxel(src.id);
        src.id = VACANT_ID;
        update_diffuseSize(src);
    }
    dest.id = pid;
    dest_mt.addVoxel(dest);
    update_diffuseSize(dest);

    return true;
}

void LatticeSpace::update_diffuseSize(Voxel& voxel)
{
    //Voxel* forward(voxel.adjoiningVoxels);
}

void LatticeSpace::remove_sparticle(ParticleID pid)
{
    if (!has_particle(pid))
    {
        return;
    }
    Voxel& voxel = voxel_as(pid);
    MolecularTypeBase* p_mt(voxel.ptr_mt);
    p_mt->removeVoxel(voxel.id);
    voxel.id = VACANT_ID;
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

Global LatticeSpace::coord2global(Integer aCoord) const
{
    Global retval;
    retval.col = aCoord / (row_size_ * layer_size_);
    retval.layer = (aCoord % (row_size_ * layer_size_)) / row_size_;
    retval.row = (aCoord % (row_size_ * layer_size_)) % row_size_;
    return retval;
}

const Position3 LatticeSpace::coord2position(Integer coord) const
{
    Integer aGlobalCol;
    Integer aGlobalLayer;
    Integer aGlobalRow;
    Global global(coord2global(coord));
    //the center point of a voxel
    Position3 position;
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            position[1] = (global.col % 2) * theHCPl + theHCPy * global.layer;
            position[2] = global.row * 2 * theNormalizedVoxelRadius +
            ((global.layer + global.col) % 2) * theNormalizedVoxelRadius;
            position[0] = global.col * theHCPx;
            break;
        case CUBIC_LATTICE:
            position[1] = global.layer * 2 * theNormalizedVoxelRadius;
            position[2] = global.row * 2 * theNormalizedVoxelRadius;
            position[0] = global.col * 2 * theNormalizedVoxelRadius;
            break;
    }
    return position;
}

/*
* Protected methods
*/
Voxel& LatticeSpace::voxel_as(ParticleID pid)
{
    for (lattice_container_type::iterator i(lattice_.begin());
            i != lattice_.end(); ++i)
    {
        if ((*i).second.id == pid)
        {
            return (*i).second;
        }
    }
    throw "Exception: Not in lattice_";
}

Voxel& LatticeSpace::voxel_at(Integer coord)
{
    lattice_container_type::iterator i(lattice_.find(coord));
    if (i == lattice_.end())
    {
        throw "Exception: Not in lattice_";
    }
    return (*i).second;
}

Integer LatticeSpace::global2coord(const Global& global) const
{
    return global.row +
        row_size_ * global.layer +
        row_size_ * layer_size_ * global.col;
}

}
