#ifndef EGFRDSIMULATOR_HPP
#define EGFRDSIMULATOR_HPP

#include "ShellID.hpp"
#include "DomainID.hpp"
#include "Shell.hpp"
#include "PairGreensFunction.hpp"
#include "ParticleSimulator.hpp"

template<typename Tworld_>
struct EGFRDSimulatorTraitsBase: public ParticleSimulatorTraitsBase<Tworld_>
{
    typedef Tworld_ world_type;
    typedef ShellID shell_id_type;
    typedef DomainID domain_id_type;
    typedef typename ParticleSimulatorTraitsBase<Tworld_>::sphere_type sphere_type;
    typedef typename ParticleSimulatorTraitsBase<Tworld_>::cylinder_type cylinder_type;
    typedef Shell<sphere_type, domain_id_type> spherical_shell_type;
    typedef Shell<cylinder_type, domain_id_type> cylindrical_shell_type;
    typedef std::pair<const shell_id_type, spherical_shell_type> spherical_shell_id_pair;
    typedef std::pair<const shell_id_type, cylindrical_shell_type> cylindrical_shell_id_pair;
    typedef int event_id_type;
    typedef EventType event_kind_type;
};

#endif /* EGFRDSIMULATOR_HPP */
