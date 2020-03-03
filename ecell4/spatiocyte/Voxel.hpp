#ifndef ECELL4_SPATIOCYTE_VOXEL_HPP
#define ECELL4_SPATIOCYTE_VOXEL_HPP

#include <boost/weak_ptr.hpp>
#include <ecell4/core/VoxelSpaceBase.hpp>
#include <ecell4/core/types.hpp>
#include <functional>

namespace ecell4
{

namespace spatiocyte
{

class SpatiocyteWorld;

struct Voxel
{
    typedef VoxelSpaceBase::coordinate_type coordinate_type;

    Voxel(boost::weak_ptr<VoxelSpaceBase> space, coordinate_type coordinate)
        : space(space), coordinate(coordinate)
    {
    }

    boost::weak_ptr<VoxelSpaceBase> space;
    coordinate_type coordinate;

public:
    const Real3 position() const
    {
        return space.lock()->coordinate2position(coordinate);
    }

    bool clear() const { return space.lock()->remove_voxel(coordinate); }

    boost::shared_ptr<VoxelPool> get_voxel_pool() const
    {
        return space.lock()->get_voxel_pool_at(coordinate);
    }

    bool operator==(const Voxel &rhs) const noexcept
    {
        return space.lock() == rhs.space.lock() && coordinate == rhs.coordinate;
    }
};

} // namespace spatiocyte

} // namespace ecell4

namespace std {
template <>
struct hash<ecell4::spatiocyte::Voxel>
{
    std::size_t operator()(const ecell4::spatiocyte::Voxel &val) const
    {
        auto ptr = val.space.lock().get();
        return hash<decltype(ptr)>()(ptr) ^
               static_cast<std::size_t>(val.coordinate);
    }
};
} // std
#endif
