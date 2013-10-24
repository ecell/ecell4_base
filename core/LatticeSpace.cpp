#include "LatticeSpace.hpp"

namespace ecell4
{

LatticeSpace::LatticeSpace()
    : VACANT_TYPE_("VACANT"), theNormalizedVoxelRadius(0.5)
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
        if (voxel.ptr_mt == &VACANT_TYPE_)
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
        Global global(coord2global(coord));
        concatenate_voxel(voxel, global);
    }
    theNullCoord = row_size_ * layer_size_ * col_size_;
}

void LatticeSpace::concatenate_voxel(Voxel& voxel, const Global& global)
{
    if (global.row > 0)
    {
        concatenate_rows(voxel, global.north());
    }
    if (global.layer > 0)
    {
        concatenate_layers(voxel, global.ventral());
    }
    if (global.col > 0)
    {
        concatenate_cols(voxel, global.west());
    }
}

void LatticeSpace::concatenate_rows(Voxel& voxel, const Global& north)
{
    Integer north_coord(global2coord(north));
    Voxel& north_voxel(voxel_at(north_coord));
    voxel.setAdjoiningVoxel(NORTH, &north_voxel);
    north_voxel.setAdjoiningVoxel(SOUTH, &voxel);
}

void LatticeSpace::concatenate_layers(Voxel& voxel, const Global& ventral)
{
    Integer ventral_coord(global2coord(ventral));
    Voxel& ventral_voxel(voxel_at(ventral_coord));
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            if ((ventral.layer + 1) % 2 + ventral.col % 2 == 1)
            {
                voxel.setAdjoiningVoxel(VENTRALN, &ventral_voxel);
                ventral_voxel.setAdjoiningVoxel(DORSALS, &voxel);
                if (ventral.row < row_size_ - 1)
                {
                    Integer ventrals_coord(global2coord(ventral.south()));
                    Voxel& ventrals_voxel(voxel_at(ventrals_coord));
                    voxel.setAdjoiningVoxel(VENTRALS, &ventrals_voxel);
                    ventrals_voxel.setAdjoiningVoxel(DORSALN, &voxel);
                }
            }
            else
            {
                voxel.setAdjoiningVoxel(VENTRALS, &ventral_voxel);
                ventral_voxel.setAdjoiningVoxel(DORSALN, &voxel);
                if (ventral.row > 0)
                {
                    Integer ventraln_coord(global2coord(ventral.north()));
                    Voxel& ventraln_voxel(voxel_at(ventraln_coord));
                    voxel.setAdjoiningVoxel(VENTRALN, &ventraln_voxel);
                    ventraln_voxel.setAdjoiningVoxel(DORSALS, &voxel);
                }
            }
            break;
        case CUBIC_LATTICE:
            voxel.setAdjoiningVoxel(VENTRAL, &ventral_voxel);
            ventral_voxel.setAdjoiningVoxel(DORSAL, &voxel);
            break;
    }
}

void LatticeSpace::concatenate_cols(Voxel& voxel, const Global& west)
{
    Integer west_coord(global2coord(west));
    Voxel& west_voxel(voxel_at(west_coord));
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            if (west.layer % 2 == 0)
            {
                if ((west.col + 1) % 2 == 1)
                {
                    voxel.setAdjoiningVoxel(NW, &west_voxel);
                    west_voxel.setAdjoiningVoxel(SE, &voxel);
                    if (west.row < row_size_ - 1)
                    {
                        Integer sw_coord(global2coord(west.south()));
                        Voxel& sw_voxel(voxel_at(sw_coord));
                        voxel.setAdjoiningVoxel(SW, &sw_voxel);
                        sw_voxel.setAdjoiningVoxel(NE, &voxel);
                    }
                    if (west.layer < layer_size_ - 1)
                    {
                        Integer c(global2coord(west.dorsal()));
                        Voxel& west_voxel(voxel_at(c));
                        voxel.setAdjoiningVoxel(WEST, &west_voxel);
                        west_voxel.setAdjoiningVoxel(EAST, &voxel);
                    }
                }
                else
                {
                    voxel.setAdjoiningVoxel(SW, &west_voxel);
                    west_voxel.setAdjoiningVoxel(NE, &voxel);
                    if (west.row > 0)
                    {
                        Integer nw_coord(global2coord(west.north()));
                        Voxel& nw_voxel(voxel_at(nw_coord));
                        voxel.setAdjoiningVoxel(NW, &nw_voxel);
                        nw_voxel.setAdjoiningVoxel(SE, &voxel);
                    }
                    if (west.layer > 0)
                    {
                        Integer c(global2coord(west.ventral()));
                        Voxel& west_voxel(voxel_at(c));
                        voxel.setAdjoiningVoxel(WEST, &west_voxel);
                        west_voxel.setAdjoiningVoxel(EAST, &voxel);
                    }
                }
            }
            else
            {
                if ((west.col + 1) % 2 == 1)
                {
                    voxel.setAdjoiningVoxel(SW, &west_voxel);
                    west_voxel.setAdjoiningVoxel(NE, &voxel);
                    if (west.row > 0)
                    {
                        Integer nw_coord(global2coord(west.north()));
                        Voxel& west_voxel(voxel_at(nw_coord));
                        voxel.setAdjoiningVoxel(NW, &west_voxel);
                        west_voxel.setAdjoiningVoxel(SE, &voxel);
                    }
                    if (west.layer < layer_size_ - 1)
                    {
                        Integer c(global2coord(west.dorsal()));
                        Voxel& west_voxel(voxel_at(c));
                        voxel.setAdjoiningVoxel(WEST, &west_voxel);
                        west_voxel.setAdjoiningVoxel(EAST, &voxel);
                    }
                }
                else
                {
                    voxel.setAdjoiningVoxel(NW, &west_voxel);
                    west_voxel.setAdjoiningVoxel(SE, &voxel);
                    if (west.row < row_size_ -1)
                    {
                        Integer sw_coord(global2coord(west.south()));
                        Voxel& sw_voxel(voxel_at(sw_coord));
                        voxel.setAdjoiningVoxel(SW, &sw_voxel);
                        sw_voxel.setAdjoiningVoxel(NE, &voxel);
                    }
                    if (west.layer > 0)
                    {
                        Integer c(global2coord(west.ventral()));
                        Voxel& west_voxel(voxel_at(c));
                        voxel.setAdjoiningVoxel(WEST, &west_voxel);
                        west_voxel.setAdjoiningVoxel(EAST, &voxel);
                    }
                }
            }
            break;
        case CUBIC_LATTICE:
            voxel.setAdjoiningVoxel(WEST, &west_voxel);
            west_voxel.setAdjoiningVoxel(EAST, &voxel);
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

    voxel0.coord = voxel1.coord;
    voxel0.diffuse_size = voxel1.diffuse_size;
    voxel0.adjoiningVoxels = voxel1.adjoiningVoxels;

    voxel1.coord = tmp_coord;
    voxel1.diffuse_size = tmp_diffuse_size;

    /*
    Global g0(coord2global(voxel0.coord));
    Global g1(coord2global(voxel1.coord));
    concatenate_voxel(voxel0, g0);
    concatenate_voxel(voxel1, g1);
    */
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
        if ((*i).second.coord == coord)
        {
            return (*i).second;
        }
    }
    throw "Exception: Not in lattice_";
}

Integer LatticeSpace::global2coord(const Global& global) const
{
    return global.row +
        row_size_ * global.layer +
        row_size_ * layer_size_ * global.col;
}

}
