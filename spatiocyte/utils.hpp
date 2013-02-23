#ifndef __ECELL4_SPATIOCYTE_UTILS_HPP
#define __ECELL4_SPATIOCYTE_UTILS_HPP

#include <ecell4/core/types.hpp>
#include <ecell4/core/Position3.hpp>

namespace ecell4
{

namespace spatiocyte
{

Real voxel_volume(const Real& voxel_radius)
{
    return 4 * std::sqrt(2) * gsl_pow_3(voxel_radius);
}

Real actual_volume_cuboid(const Position3& edge_lengths, const Real& voxel_radius)
{
    const Real normalized_voxel_radius(0.5);
    // const Integer hcpl(normalized_voxel_radius / std::sqrt(3));
    const Real hcpx(normalized_voxel_radius * std::sqrt(8.0 / 3.0)),
        hcpy(normalized_voxel_radius * std::sqrt(3));
    const Position3 center(divide(edge_lengths, 2 * (2 * voxel_radius)));
    Integer num_rows(
        static_cast<Integer>(rint(center[2] / normalized_voxel_radius)) + 2),
        num_layers(static_cast<Integer>(rint(center[1] * 2 / hcpy)) + 2),
        num_cols(static_cast<Integer>(rint(center[0] * 2 / hcpx)) + 2);

    // readjust_surface_boundary_sizes
    num_cols += (num_cols % 2 == 0) ? 0 : 1;
    num_layers += (num_layers % 2 == 0) ? 0 : 1;
    num_rows += (num_rows % 2 == 0) ? 0 : 1;

    if (false) // isPeriodicEdge
    {
        num_cols += 1;
        num_layers += 1;
        num_rows += 1;
    }

    const Integer num_voxels(num_rows * num_layers * num_cols);
    return voxel_volume(voxel_radius) * num_voxels;
}

} // spatiocyte

} // ecell4

#endif /* __ECELL4_SPATIOCYTE_UTILS_HPP */
