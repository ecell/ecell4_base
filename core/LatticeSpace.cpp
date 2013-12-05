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

LatticeSpace::LatticeSpace(const Position3& edge_lengths)
{
    this->edge_lengths_ = edge_lengths;
    set_lattice_properties();
}

LatticeSpace::~LatticeSpace()
{
}

Integer LatticeSpace::num_species() const
{
    return spmap_.size();
}

bool LatticeSpace::has_species(const Species& sp) const
{
    return spmap_.find(sp) != spmap_.end();
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
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        count += (*itr).second.voxels().size();
    }

    return count;
}

Integer LatticeSpace::num_particles(const Species& sp) const
{
    spmap::const_iterator itr(spmap_.find(sp));
    Integer count(0);
    if (itr != spmap_.end())
    {
        count = (*itr).second.voxels().size();
    }
    return count;
}

bool LatticeSpace::has_particle(const ParticleID& pid) const
{
    bool flg(false);
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularType& mt((*itr).second);
        for (MolecularType::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            if ((*vitr).second == pid)
            {
                flg = true;
                break;
            }
        }
    }

    return flg;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles() const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularType& mt((*itr).second);
        for (MolecularType::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            retval.push_back(std::pair<ParticleID, Particle>(
                        (*vitr).second, particle_at((*vitr).first)));
        }
    }

    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    spmap::const_iterator itr(spmap_.find(sp));
    if (itr != spmap_.end())
    {
        const MolecularType& mt((*itr).second);
        for (MolecularType::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            retval.push_back(std::pair<ParticleID, Particle>(
                        (*vitr).second, particle_at((*vitr).first)));
        }
    }

    return retval;
}

std::vector<SParticle> LatticeSpace::list_sparticles() const
{
    std::vector<SParticle> retval;
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularType& mt((*itr).second);
        std::vector<SParticle> sparticles(mt.sparticles());
        retval.insert(retval.end(), sparticles.begin(), sparticles.end());
    }

    return retval;
}

std::vector<SParticle> LatticeSpace::list_sparticles(const Species& sp) const
{
    spmap::const_iterator itr(spmap_.find(sp));
    if (itr != spmap_.end())
    {
        const MolecularType& mt((*itr).second);
        return mt.sparticles();
    }

    std::vector<SParticle> retval;
    return retval;
}

bool LatticeSpace::update_particle(const ParticleID& pid, const Particle& p)
{
    Integer coord(position2coord(p.position()));
    MolecularType* dest_mt = get_molecular_type(p.species());

    if (has_particle(pid))
    {
        Integer coord2(get_coord(pid));
        MolecularTypeBase* src_ptr_mt(voxels_.at(coord));
        src_ptr_mt->removeVoxel(coord2);
        voxel_container::iterator itr(voxels_.begin() + coord2);
        voxels_.erase(itr);
        voxels_.insert(itr, NULL);
    }
    dest_mt->addVoxel(coord, pid);

    voxel_container::iterator itr(voxels_.begin() + coord);
    voxels_.erase(itr);
    voxels_.insert(itr, dest_mt);

    return true;
}

bool LatticeSpace::move(Integer from, Integer to)
{
    if (has_particle(to) || has_particle(from))
        return false;

    MolecularType* mt(voxel_at(from));
    MolecularType::container_type::iterator mitr(mt->find(from));
    if (mitr == mt->end())
        return false;

    ParticleID pid((*mitr).second);
    mt->removeVoxel(from);
    mt->addVoxel(to, pid);

    voxel_container::iterator itr;
    itr = voxels_.begin() + from;
    voxels_.erase(itr);
    voxels_.insert(itr, NULL);
    itr = voxels_.begin() + to;
    voxels_.erase(itr);
    voxels_.insert(itr, mt);

    return true;
}

bool LatticeSpace::has_particle(Integer coord) const
{
    return voxels_.at(coord) != NULL;
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

    Real lengthX = edge_lengths_[0];
    Real lengthY = edge_lengths_[1];
    Real lengthZ = edge_lengths_[2];

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
        voxels_.push_back(NULL);
    }
}

/*
 * original methods
 */

MolecularType* LatticeSpace::get_molecular_type(const Species& sp)
{
    spmap::iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        MolecularType mt(sp);
        spmap_.insert(spmap::value_type(sp, mt));
        itr = spmap_.find(sp);
    }
    return &((*itr).second);
}

const Global LatticeSpace::position2global(const Position3& pos) const
{
    Global global;
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            global.col = (Integer)(pos[0] / theHCPx);
            global.layer = (Integer)((pos[1] - (global.col % 2) * theHCPl) / theHCPy);
            global.row = (Integer)((pos[2] - (((global.layer + global.col) % 2) *
                            theNormalizedVoxelRadius)) /
                    theNormalizedVoxelRadius / 2);
            break;

        case CUBIC_LATTICE:
            global.layer = (Integer)(pos[1] / theNormalizedVoxelRadius / 2);
            global.row = (Integer)(pos[2] / theNormalizedVoxelRadius / 2);
            global.col = (Integer)(pos[0] / theNormalizedVoxelRadius / 2);
            break;
    }
    return global;
}

const Global LatticeSpace::coord2global(Integer aCoord) const
{
    Global retval;
    retval.col = aCoord / (row_size_ * layer_size_);
    retval.layer = (aCoord % (row_size_ * layer_size_)) / row_size_;
    retval.row = (aCoord % (row_size_ * layer_size_)) % row_size_;
    return retval;
}

const Position3 LatticeSpace::global2position(const Global& global) const
{
    return coord2position(global2coord(global));
}

const Position3 LatticeSpace::coord2position(Integer coord) const
{
    const Global global(coord2global(coord));
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

Integer LatticeSpace::position2coord(const Position3& pos) const
{
    return global2coord(position2global(pos));
}

Integer LatticeSpace::get_coord(const ParticleID& pid) const
{
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularType& mt((*itr).second);
        for (MolecularType::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            if ((*vitr).second == pid)
            {
                return (*vitr).first;
            }
        }
    }
    throw "Exception: Not in lattice";
}

MolecularType* LatticeSpace::voxel_at(Integer coord) const
{
    return voxels_.at(coord);
}

Integer LatticeSpace::global2coord(const Global& global) const
{
    return global.row +
        row_size_ * global.layer +
        row_size_ * layer_size_ * global.col;
}

const Particle LatticeSpace::particle_at(Integer coord) const
{
    const MolecularType* ptr_mt(voxel_at(coord));
    const Species& sp = ptr_mt->species();
    const Position3& pos = coord2position(coord);
    const Real& radius = 0;
    const Real& D = 0;
    Particle particle(sp, pos, radius, D);
    return particle;
}

} // ecell4
