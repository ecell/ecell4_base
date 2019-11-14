#include "HCPLatticeSpace.hpp"

#include <cmath>

#ifdef WIN32_MSC
#include <boost/numeric/interval/detail/msvc_rounding_control.hpp>
#endif

namespace ecell4
{

#ifdef WIN32_MSC
double rint(const double x)
{
    return boost::numeric::interval_lib::detail::rint(x);
}

double round(const double x) { return floor(x + 0.5); }
#endif

void HCPLatticeSpace::set_lattice_properties(const Real3 &edge_lengths,
                                             const bool is_periodic)
{
    // XXX: derived from SpatiocyteStepper::setLatticeProperties()
    HCP_L = voxel_radius_ / sqrt(3.0);
    HCP_X = voxel_radius_ * sqrt(8.0 / 3.0); // Lx
    HCP_Y = voxel_radius_ * sqrt(3.0);       // Ly

    const Real &lengthX = edge_lengths[0];
    const Real &lengthY = edge_lengths[1];
    const Real &lengthZ = edge_lengths[2];

    col_size_ = (Integer)rint(lengthX / HCP_X);
    layer_size_ = (Integer)rint(lengthY / HCP_Y);
    row_size_ = (Integer)rint((lengthZ / 2) / voxel_radius_);

    if (is_periodic)
    {
        // The number of voxels in each axis must be even for a periodic
        // boundary.
        col_size_ = (col_size_ % 2 == 0 ? col_size_ : col_size_ + 1);
        layer_size_ = (layer_size_ % 2 == 0 ? layer_size_ : layer_size_ + 1);
        row_size_ = (row_size_ % 2 == 0 ? row_size_ : row_size_ + 1);
    }

    edge_lengths_ = Real3(col_size_ * HCP_X, layer_size_ * HCP_Y,
                          row_size_ * voxel_radius_ * 2);

    row_size_ += 2;
    layer_size_ += 2;
    col_size_ += 2;
}

} // namespace ecell4
