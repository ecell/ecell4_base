#ifndef TRAITS_HPP
#define TRAITS_HPP

#include <exception>
#include <stdexcept>

//#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <boost/range/size.hpp>

#include "../utils/array_traits.hpp"
#include "../utils.hpp"
#include "../geometry.hpp"
#include "../MatrixSpace.hpp"
#include "../Vector3.hpp"
#include "../Sphere.hpp"
#include "../Cylinder.hpp"
#include "../Box.hpp"
#include "../Plane.hpp"
#include "../Point.hpp"
#include "../Model.hpp"
#include "../World.hpp"
#include "../Multi.hpp"
#include "../ShapedDomain.hpp"
#include "../GSLRandomNumberGenerator.hpp"
#include "../EGFRDSimulator.hpp"
#include "../BDPropagator.hpp"
#include "../StructureUtils.hpp"
#include "../AnalyticalSingle.hpp"
#include "../AnalyticalPair.hpp"
#include "../EventScheduler.hpp"

#include "peer/wrappers/range/pyiterable_range.hpp"

namespace binding {

typedef ::not_found NotFound;
typedef ::already_exists AlreadyExists;
typedef ::illegal_state IllegalState;
typedef ::GSLRandomNumberGenerator GSLRandomNumberGenerator;
typedef ::CyclicWorldTraits<Real, Real> WorldTraits;
typedef WorldTraits::particle_type Particle;
typedef WorldTraits::structure_id_type StructureID;
typedef WorldTraits::species_id_type SpeciesID;
typedef WorldTraits::species_type SpeciesInfo;
typedef WorldTraits::structure_type Structure;
typedef WorldTraits::length_type Length;
typedef WorldTraits::position_type Position;
typedef WorldTraits::surface_type Surface;
typedef WorldTraits::region_type Region;
typedef ::World<WorldTraits> World;
typedef ::Model Model;
typedef ::NetworkRules NetworkRules;
typedef NetworkRules::reaction_rule_generator ReactionRuleGenerator;
typedef World::transaction_type Transaction;
typedef World::base_type::base_type ParticleContainer; 
typedef ::TransactionImpl<ParticleContainer> TransactionImpl;
typedef ::EGFRDSimulatorTraitsBase<World> EGFRDSimulatorTraits;
typedef ::EGFRDSimulator<EGFRDSimulatorTraits> EGFRDSimulator;
typedef ::MultiParticleContainer<EGFRDSimulatorTraits> MultiParticleContainer;
typedef EGFRDSimulatorTraits::box_type Box;
typedef EGFRDSimulatorTraits::sphere_type Sphere;
typedef EGFRDSimulatorTraits::cylinder_type Cylinder;
typedef EGFRDSimulatorTraits::plane_type Plane;
typedef ::BDPropagator<EGFRDSimulatorTraits> BDPropagator;
typedef EGFRDSimulatorTraits::shell_id_type ShellID;
typedef EGFRDSimulatorTraits::domain_id_type DomainID;
typedef EGFRDSimulatorTraits::reaction_rule_type ReactionRuleInfo;
typedef EGFRDSimulatorTraits::network_rules_type NetworkRulesWrapper;
typedef EGFRDSimulator::domain_type Domain;
typedef ::ShapedDomain<EGFRDSimulatorTraits> ShapedDomain;
typedef EGFRDSimulator::spherical_shell_type SphericalShell;
typedef EGFRDSimulator::cylindrical_shell_type CylindricalShell;
typedef EGFRDSimulator::spherical_single_type::base_type Single;
typedef EGFRDSimulator::spherical_pair_type::base_type Pair;
typedef EGFRDSimulator::spherical_single_type SphericalSingle;
typedef EGFRDSimulator::cylindrical_single_type CylindricalSingle;
typedef EGFRDSimulator::spherical_pair_type SphericalPair;
typedef EGFRDSimulator::cylindrical_pair_type CylindricalPair;
typedef EGFRDSimulator::multi_type Multi;
typedef ::MatrixSpace<SphericalShell, ShellID> SphericalShellContainer;
typedef ::MatrixSpace<CylindricalShell, ShellID> CylindricalShellContainer;
typedef ::StructureUtils<EGFRDSimulatorTraits> StructureUtils;
typedef EGFRDSimulatorTraits::planar_surface_type PlanarSurface;
typedef EGFRDSimulatorTraits::spherical_surface_type SphericalSurface;
typedef EGFRDSimulatorTraits::cylindrical_surface_type CylindricalSurface;
typedef EGFRDSimulatorTraits::cuboidal_region_type CuboidalRegion;
typedef EGFRDSimulatorTraits::reaction_record_type ReactionRecord;
typedef EGFRDSimulatorTraits::reaction_recorder_type ReactionRecorder;
} // namespace binding

#endif /* TRAITS_HPP */
