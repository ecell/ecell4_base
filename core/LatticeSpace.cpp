#include "LatticeSpace.hpp"

namespace ecell4
{

LatticeSpace::LatticeSpace()
    : VACANT_TYPE_("VACANT")
{
    set_lattice_properties();
    construct_lattice();
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
    return lattice_.size();
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

// TODO
void LatticeSpace::set_lattice_properties()
{
//    stride_ = UINT_MAX / get_species_set().size();
    stride_ = 100000000 / get_species_set().size();
//    Comp* rootComp(theComps[0]);
    Comp* rootComp;
//    rotateCompartment(rootComp);
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            adjoining_size_ = 12;
            theHCPl = theNormalizedVoxelRadius/sqrt(3); 
            theHCPx = theNormalizedVoxelRadius*sqrt(8.0/3); //Lx
            theHCPy = theNormalizedVoxelRadius*sqrt(3); //Ly
            break;
        case CUBIC_LATTICE:
            adjoining_size_ = 6;
            break;
    }
    if(rootComp->geometry == CUBOID)
    {
        //We do not give any leeway space between the simulation boundary
        //and the cell boundary if it is CUBOID to support
        //periodic boundary conditions:
        theCenterPoint[2] = rootComp->lengthZ / 2; //row
        theCenterPoint[1] = rootComp->lengthY / 2; //layer
        theCenterPoint[0] = rootComp->lengthX / 2; //column
    }
    else
    {
        switch(lattice_type_)
        {
            case HCP_LATTICE:
                theCenterPoint[2] = rootComp->lengthZ / 2 + 4 *
                theNormalizedVoxelRadius; //row
                theCenterPoint[1] = rootComp->lengthY / 2 + 2 * theHCPy; //layer
                theCenterPoint[0] = rootComp->lengthX / 2 + 2 * theHCPx; //column
                break;
            case CUBIC_LATTICE:
                theCenterPoint[2] = rootComp->lengthZ / 2 +
                    8 * theNormalizedVoxelRadius;
                theCenterPoint[1] = rootComp->lengthY / 2 +
                    8 * theNormalizedVoxelRadius;
                theCenterPoint[0] = rootComp->lengthX / 2 +
                    8 * theNormalizedVoxelRadius;
                break;
        }
    }
    rootComp->centerPoint = theCenterPoint;
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            row_size_ = (Integer)rint((theCenterPoint[2])/
                                      (theNormalizedVoxelRadius));
            layer_size_ = (Integer)rint((theCenterPoint[1]*2)/theHCPy);
            col_size_ = (Integer)rint((theCenterPoint[0]*2)/theHCPx);
            break;
        case CUBIC_LATTICE:
            row_size_ = (Integer)rint((theCenterPoint[2])/
                                      (theNormalizedVoxelRadius));
            layer_size_ = (Integer)rint((theCenterPoint[1])/
                                      (theNormalizedVoxelRadius));
            col_size_ = (Integer)rint((theCenterPoint[0])/
                                      (theNormalizedVoxelRadius));
            break;
    }
    //For the CUBOID cell geometry, we need to readjust the size of
    //row, layer and column according to the boundary condition of its surfaces
    //to reflect the correct volume. This is because periodic boundary will
    //consume a layer of the surface voxels:
    if(rootComp->geometry == CUBOID)
    {
        //We need to increase the row, layer and col size by 2 because
        //the entire volume must be surrounded by nullID voxels to avoid
        //self-homodimerization reaction.
        row_size_ += 2;
        col_size_ += 2;
        layer_size_ += 2;
//        readjustSurfaceBoundarySizes();
    }

    //You should only resize lattice_ once here. You should never
    //push or readjust the capacity of lattice_ since all the pointers
    //to the voxel will become invalid once you do that:
    //Also add one more voxel for the nullVoxel:
    lattice_.resize(row_size_ * layer_size_ * col_size_ + 1);
    //Initialize the null coord:
    theNullCoord = row_size_ * layer_size_ * col_size_;
    //lattice_[theNullCoord].idx = theNullID * stride_;
}

// TODO
void LatticeSpace::construct_lattice()
{
//    Comp* rootComp(theComps[0]);
    Comp* rootComp;
    Integer size(row_size_ * layer_size_ * col_size_);
    Integer a = 0;
//    unsigned short rootID(rootComp->vacantSpecies->getID());
    Integer rootID = 1; // TODO
    for (lattice_container_type::iterator i(lattice_.begin());
            a != size; ++i, ++a)
    {
        (*i).coord = a;
        (*i).adjoiningVoxels.resize(adjoining_size_);
        (*i).diffuse_size = adjoining_size_;
        Integer col = a / (row_size_ * layer_size_);
        Integer layer = (a % (row_size_ * layer_size_)) / row_size_;
        Integer row = (a % (row_size_ * layer_size_)) % row_size_;
//        if(rootComp->geometry == CUBOID || is_inside_coord(a, rootComp, 0))
        if(is_inside_coord(a, rootComp, 0))
        {
            //By default, the voxel is vacant and we set it to the root id:
//            (*i).idx = rootID*stride_;
            for(Integer j = 0; j != adjoining_size_; ++j)
            {
                // By default let the adjoining voxel pointer point to the
                // source voxel (i.e., itself)
                (*i).adjoiningVoxels[j] = &*i;
            }
            concatenate_voxel(*i, row, layer, col);
        }
        else
        {
            //We set id = theNullID if it is an invalid voxel, i.e., no molecules
            //will occupy it:
//            (*i).idx = theNullID*stride_;
            //Concatenate some of the null voxels close to the surface:
//            if(is_inside_coord(a, rootComp, 4))
            if (is_inside_coord(a, rootComp, 4))
            {
                concatenate_voxel(*i, row, layer, col);
            }
        }
    }
//    if(rootComp->geometry == CUBOID)
//    {
//        concatenatePeriodicSurfaces();
//    }
}

void LatticeSpace::concatenate_voxel(Voxel& voxel,
        Integer row, Integer layer, Integer col)
{
    Integer coord = position2coord(row, layer, col);
    if (row > 0)
    {
        concatenate_rows(voxel, row-1, layer, col);
    }
    if (layer > 0)
    {
        concatenate_layers(voxel, row, layer-1, col);
    }
    if (col > 0)
    {
        concatenate_cols(voxel, row, layer, col-1);
    }
}

void LatticeSpace::concatenate_rows(Voxel& voxel,
        Integer row, Integer layer, Integer col)
{
    Integer b(position2coord(row, layer, col));
    Voxel& adjoining(lattice_.at(b));
    voxel.adjoiningVoxels[NORTH] = &adjoining;
    adjoining.adjoiningVoxels[SOUTH] = &voxel;
}

void LatticeSpace::concatenate_layers(Voxel& voxel,
        Integer row, Integer layer, Integer col)
{
    Integer b(position2coord(row, layer, col));
    Voxel& adjoining(lattice_.at(b));
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            if((layer + 1) % 2 + col % 2 == 1)
            {
                voxel.adjoiningVoxels[VENTRALN] = &adjoining;
                adjoining.adjoiningVoxels[DORSALS] = &voxel;
                if(row < row_size_-1)
                {
                    Integer c(position2coord(row + 1, layer, col));
                    Voxel& adjoining(lattice_.at(c));
                    voxel.adjoiningVoxels[VENTRALS] = &adjoining;
                    adjoining.adjoiningVoxels[DORSALN] = &voxel;
                }
            }
            else
            {
                voxel.adjoiningVoxels[VENTRALS] = &adjoining;
                adjoining.adjoiningVoxels[DORSALN] = &voxel;
                if(row > 0)
                {
                    Integer c(position2coord(row - 1, layer, col));
                    Voxel& adjoining(lattice_.at(c));
                    voxel.adjoiningVoxels[VENTRALN] = &adjoining;
                    adjoining.adjoiningVoxels[DORSALS] = &voxel;
                }
            }
            break;
        case CUBIC_LATTICE:
            voxel.adjoiningVoxels[VENTRAL] = &adjoining;
            adjoining.adjoiningVoxels[DORSAL] = &voxel;
            break;
    }
}

void LatticeSpace::concatenate_cols(Voxel& voxel,
        Integer row, Integer layer, Integer col)
{
    Integer b(position2coord(row, layer, col));
    Voxel& adjoining(lattice_.at(b));
    switch(lattice_type_)
    {
        case HCP_LATTICE:
            if(layer % 2 == 0)
            {
                if((col + 1) % 2 == 1)
                {
                    voxel.adjoiningVoxels[NW] = &adjoining;
                    adjoining.adjoiningVoxels[SE] = &voxel;
                    if(row < row_size_ - 1)
                    {
                        Integer c(position2coord(row + 1, layer, col));
                        Voxel& adjoining(lattice_.at(c));
                        voxel.adjoiningVoxels[SW] = &adjoining;
                        adjoining.adjoiningVoxels[NE] = &voxel;
                    }
                    if(layer < layer_size_-1)
                    {
                        Integer c(position2coord(row, layer + 1, col));
                        Voxel& adjoining(lattice_.at(c));
                        voxel.adjoiningVoxels[WEST] = &adjoining;
                        adjoining.adjoiningVoxels[EAST] = &voxel;
                    }
                }
                else
                {
                    voxel.adjoiningVoxels[SW] = &adjoining;
                    adjoining.adjoiningVoxels[NE] = &voxel;
                    if(row > 0)
                    {
                        Integer c(position2coord(row - 1, layer, col));
                        Voxel& adjoining(lattice_.at(c));
                        voxel.adjoiningVoxels[NW] = &adjoining;
                        adjoining.adjoiningVoxels[SE] = &voxel;
                    }
                    if(layer > 0)
                    {
                        Integer c(position2coord(row, layer - 1, col));
                        Voxel& adjoining(lattice_.at(c));
                        voxel.adjoiningVoxels[WEST] = &adjoining;
                        adjoining.adjoiningVoxels[EAST] = &voxel;
                    }
                }
            }
            else
            {
                if((col + 1) % 2 == 1)
                {
                    voxel.adjoiningVoxels[SW] = &adjoining;
                    adjoining.adjoiningVoxels[NE] = &voxel;
                    if(row > 0)
                    {
                        Integer c(position2coord(row - 1, layer, col));
                        Voxel& adjoining(lattice_.at(c));
                        voxel.adjoiningVoxels[NW] = &adjoining;
                        adjoining.adjoiningVoxels[SE] = &voxel;
                    }
                    if(layer < layer_size_-1)
                    {
                        Integer c(position2coord(row, layer + 1, col));
                        Voxel& adjoining(lattice_.at(c));
                        voxel.adjoiningVoxels[WEST] = &adjoining;
                        adjoining.adjoiningVoxels[EAST] = &voxel;
                    }
                }
                else
                {
                    voxel.adjoiningVoxels[NW] = &adjoining;
                    adjoining.adjoiningVoxels[SE] = &voxel;
                    if(row < row_size_ - 1)
                    {
                        Integer c(position2coord(row + 1, layer, col));
                        Voxel& adjoining(lattice_.at(c));
                        voxel.adjoiningVoxels[SW] = &adjoining;
                        adjoining.adjoiningVoxels[NE] = &voxel;
                    }
                    if(layer > 0)
                    {
                        Integer c(position2coord(row, layer - 1, col));
                        Voxel& adjoining(lattice_.at(c));
                        voxel.adjoiningVoxels[WEST] = &adjoining;
                        adjoining.adjoiningVoxels[EAST] = &voxel;
                    }
                }
            }
            break;
        case CUBIC_LATTICE:
            voxel.adjoiningVoxels[WEST] = &adjoining;
            adjoining.adjoiningVoxels[EAST] = &voxel;
            break;
    }
}

// TODO
bool LatticeSpace::is_inside_coord(Integer coord, Comp* comp,
        Real delta)
{
/*
    Point aPoint(coord2point(aCoord));
    Point aCenterPoint(aComp->centerPoint);
    Point aWestPoint(aComp->centerPoint);
    Point anEastPoint(aComp->centerPoint);
    aPoint.x -= aCenterPoint.x;
    aPoint.y -= aCenterPoint.y;
    aPoint.z -= aCenterPoint.z;
    rotateX(aComp->rotateX, &aPoint);
    rotateY(aComp->rotateY, &aPoint);
    rotateZ(aComp->rotateZ, &aPoint);
    aPoint.x += aCenterPoint.x;
    aPoint.y += aCenterPoint.y;
    aPoint.z += aCenterPoint.z;
    double aRadius(0);
    switch(aComp->geometry)
    {
    case CUBOID:
        if(sqrt(pow(aPoint.x-aCenterPoint.x, 2)) <=
        aComp->lengthX/2+theNormalizedVoxelRadius+delta &&
        sqrt(pow(aPoint.y-aCenterPoint.y, 2)) <=
        aComp->lengthY/2+theNormalizedVoxelRadius+delta &&
        sqrt(pow(aPoint.z-aCenterPoint.z, 2)) <=
        aComp->lengthZ/2+theNormalizedVoxelRadius+delta)
        {
            return true;
        }
        break;
    case ELLIPSOID:
        //If the distance between the voxel and the center point is less than
        //or equal to radius-2, then the voxel cannot be a surface voxel:
        if(pow(aPoint.x-aCenterPoint.x, 2)/pow((aComp->lengthX+delta)/2, 2)+
        pow(aPoint.y-aCenterPoint.y, 2)/pow((aComp->lengthY+delta)/2, 2)+
        pow(aPoint.z-aCenterPoint.z, 2)/pow((aComp->lengthZ+delta)/2, 2) <= 1)
        {
            return true;
        }
        break;
    case CYLINDER:
        //The axial point of the cylindrical portion of the rod:
        aCenterPoint.x = aPoint.x;
        aWestPoint.x = aComp->centerPoint.x-aComp->lengthX/2;
        anEastPoint.x = aComp->centerPoint.x+aComp->lengthX/2;
        aRadius = aComp->lengthY/2+theNormalizedVoxelRadius;
        //If the distance between the voxel and the center point is less than
        //or equal to the radius, then the voxel must be inside the Comp:
        if((aPoint.x >= aWestPoint.x && aPoint.x <= anEastPoint.x &&
        distance(aPoint, aCenterPoint) <= aRadius+delta))
        {
            return true;
        }
        break;
    case ROD:
        //The axial point of the cylindrical portion of the rod:
        aCenterPoint.x = aPoint.x;
        aWestPoint.x = aComp->centerPoint.x-aComp->lengthX/2+aComp->lengthY/2;
        anEastPoint.x = aComp->centerPoint.x+aComp->lengthX/2-aComp->lengthY/2;
        aRadius = aComp->lengthY/2+theNormalizedVoxelRadius;
        //If the distance between the voxel and the center point is less than
        //or equal to the radius, then the voxel must be inside the Comp:
        if((aPoint.x >= aWestPoint.x && aPoint.x <= anEastPoint.x &&
        distance(aPoint, aCenterPoint) <= aRadius+delta) ||
        (aPoint.x < aWestPoint.x &&
        distance(aPoint, aWestPoint) <= aRadius+delta) ||
        (aPoint.x > anEastPoint.x &&
        distance(aPoint, anEastPoint) <= aRadius+delta))
        {
            return true;
        }
        break;
    case PYRAMID:
        aRadius = ((aCenterPoint.y+aComp->lengthY/2)-aPoint.y)/aComp->lengthY;
        if(sqrt(pow(aPoint.y-aCenterPoint.y, 2)) <= aComp->lengthY/2+delta &&
        sqrt(pow(aPoint.x-aCenterPoint.x, 2)) <=
        aComp->lengthX*aRadius/2+delta &&
        sqrt(pow(aPoint.z-aCenterPoint.z, 2)) <=
        aComp->lengthZ*aRadius/2+delta)
        {
            return true;
        }
        break;
    case ERYTHROCYTE:
        if(delta > 0)
        {
            return true;
        }
        else if(delta < 0)
        {
            return false;
        }
        const double Rsq(pow(aPoint.x-aCenterPoint.x, 2)/
        pow((aComp->lengthX)/2, 2)+
        pow(aPoint.y-aCenterPoint.y, 2)/
        pow((aComp->lengthY)/2, 2));
        if(Rsq > 1)
        {
            return false;
        }
        const double a(0.5);
        const double b(0.1);
        const double R(sqrt(Rsq));
        const double thickness(((1-cos(M_PI*0.5*R))*(a-b)+b)*sqrt(1-Rsq));
        const double height((aPoint.z-aCenterPoint.z)/(2*(aComp->lengthZ)));
        if(thickness*thickness >= height*height)
        {
            return true;
        }
        break;
    }
*/
    return false;
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
    dest.p_molecule_type = src.p_molecule_type;
    dest.id = src.id;
    src.p_molecule_type = &VACANT_TYPE_;
    src.id = VACANT_ID_;
    return true;
}

void LatticeSpace::remove_sparticle(ParticleID pid)
{
    Voxel& v = voxel_as(pid);
    v.p_molecule_type = &VACANT_TYPE_;
    v.id = VACANT_ID_;
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
        Species s = (*i).species();
        retval.insert(&s);
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
