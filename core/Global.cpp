#include "Global.hpp"

namespace ecell4 {

Global Global::east() const
{
    Global retval(*this);
    retval.col += 1;
    return retval;
}

Global Global::west() const
{
    Global retval(*this);
    retval.col -= 1;
    return retval;
}

Global Global::south() const
{
    Global retval(*this);
    retval.row += 1;
    return retval;
}

Global Global::north() const
{
    Global retval(*this);
    retval.row -= 1;
    return retval;
}

Global Global::dorsal() const
{
    Global retval(*this);
    retval.layer += 1;
    return retval;
}

Global Global::ventral() const
{
    Global retval(*this);
    retval.layer -= 1;
    return retval;
}

} // ecell4
