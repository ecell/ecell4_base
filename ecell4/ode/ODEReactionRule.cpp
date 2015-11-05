
#include "ODEReactionRule.hpp"
#include <sstream>
#include <boost/format.hpp>

namespace ecell4
{

namespace ode
{

const std::string ODEReactionRule::as_string() const
{
    std::stringstream ss_left_side, ss_right_side, ss_k_side;

    bool first = true;
    for (reaction_leftside_container_type::const_iterator it(reactants_.begin());
         it != reactants_.end(); it++)
    {
        if (!first)
        {
            ss_left_side << "+";
        }
        else
        {
            first = false;
        }

        if (it->first == 1)
        {
            ss_left_side << boost::format("%s") % it->second.serial();
        }
        else
        {
            ss_left_side << boost::format("%g*%s") % it->first % it->second.serial();
        }
    }

    first = true;
    for (reaction_rightside_container_type::const_iterator it(products_.begin());
         it != products_.end(); it++)
    {
        if (!first)
        {
            ss_right_side << "+";
        }
        else
        {
            first = false;
        }

        if (it->first == 1)
        {
            ss_right_side << boost::format("%s") % it->second.serial();
        }
        else
        {
            ss_right_side << boost::format("%g*%s") % it->first % it->second.serial();
        }
    }

    if (!this->has_ratelaw())
    {
        ss_k_side << "nan";
    }
    else
    {
        ss_k_side << this->ratelaw_->as_string();
    }
    return (boost::format("%s>%s|%g") % ss_left_side.str() % ss_right_side.str() % ss_k_side.str()).str();
}

}   // ode

}   // ecell4
