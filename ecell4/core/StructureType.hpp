#ifndef ECELL4_STRUCTURE_TYPE_HPP
#define ECELL4_STRUCTURE_TYPE_HPP

#include "VoxelPool.hpp"


namespace ecell4
{

class StructureType
    : public VoxelPool
{
private:

    typedef VoxelPool base_type;

public:

    StructureType(
        const Species& species, boost::weak_ptr<VoxelPool> location,
        const Real& radius = 0.0)
        : base_type(species, location, radius, 0),
          size_(0)
    {
        ;
    }

    virtual ~StructureType()
    {
        ;
    }

    virtual voxel_type_type const voxel_type() const
    {
        return STRUCTURE;
    }

    const Integer size() const
    {
        return size_;
    }

    void add_voxel(const coordinate_id_pair_type& info)
    {
        if (info.pid != ParticleID())
        {
            throw NotSupported("No ParticleID is allowed.");
        }

        ++size_;
    }

    coordinate_id_pair_type pop(const coordinate_type& coord)
    {
        --size_;
        return coordinate_id_pair_type(ParticleID(), coord);
    }

    bool remove_voxel_if_exists(const coordinate_type& coord)
    {
        --size_;
        return true;
    }

private:

    Integer size_;

};

} //ecell4

#endif /* ECELL4_STRUCTURE_TYPE_HPP */
