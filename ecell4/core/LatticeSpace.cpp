#include "LatticeSpace.hpp"
// #include <cmath>
// #include <sstream>
// #include <algorithm>

namespace ecell4
{

bool LatticeSpace::can_move(const coordinate_type& src,
        const coordinate_type& dest) const
{
    return false;
}

bool LatticeSpace::make_structure_type(const Species& sp,
    Shape::dimension_kind dimension, const std::string loc)
{
    return false;
}

bool LatticeSpace::make_interface_type(const Species& sp,
    Shape::dimension_kind dimension, const std::string loc)
{
    return false;
}

} // ecell4
