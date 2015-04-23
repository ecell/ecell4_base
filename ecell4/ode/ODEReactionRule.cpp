
#include "ODEReactionRule.hpp"
#include <sstream>
#include <boost/format.hpp>

namespace ecell4 
{
namespace ode
{

const std::string ODEReactionRule::as_string() const
{
    std::stringstream ss_left_side, ss_right_side;
    bool first = true;
    for(reaction_leftside_container_type::const_iterator it(reactants_.begin());
            it != reactants_.end(); it++)
    {
        if (first == false)
        {
            ss_left_side << boost::format("+ %f(%s) ") % it->first % it->second.serial();
        }
        else
        {
            ss_left_side << boost::format("%f(%s) ") % it->first % it->second.serial();
            first = false;
        }
    }
    first = true;
    for(reaction_rightside_container_type::const_iterator it(products_.begin());
            it != products_.end(); it++)
    {
        if (first == false)
        {
            ss_right_side << boost::format("+ %f(%s) ") % it->first % it->second.serial();
        }
        else
        {
            ss_right_side << boost::format("%f(%s) ") % it->first % it->second.serial();
            first = false;
        }
    }
    return (boost::format("%s ---> %s  (k: %f)") % ss_left_side.str() % ss_right_side.str() % k() ).str();
}

}   // ode

}   // ecell4
