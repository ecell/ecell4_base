#ifndef EGFRDSIMULATOR_HPP
#define EGFRDSIMULATOR_HPP

#include "ShellID.hpp"
#include "DomainID.hpp"
#include "SphericalShell.hpp"
#include "CylindricalShell.hpp"

template<typename Tworld_>
struct EGFRDSimulatorTraitsBase
{
    typedef Tworld_ world_type;
    typedef ShellID shell_id_type;
    typedef DomainID domain_id_type;
    typedef SphericalShell<typename world_type::length_type,
                           domain_id_type> spherical_shell_type;
    typedef CylindricalShell<typename world_type::length_type,
                             domain_id_type> cylindrical_shell_type;
};

#endif /* EGFRDSIMULATOR_HPP */
