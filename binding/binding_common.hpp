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
#include "../ParticleModel.hpp"
#include "../World.hpp"
#include "../Multi.hpp"
#include "../ShapedDomain.hpp"
#include "../GSLRandomNumberGenerator.hpp"
#include "../ParticleSimulator.hpp"
#include "../EGFRDSimulator.hpp"
#include "../BDPropagator.hpp"
#include "../BDSimulator.hpp"
#include "../StructureUtils.hpp"
#include "../AnalyticalSingle.hpp"
#include "../AnalyticalPair.hpp"
#include "../EventScheduler.hpp"
#include "../Logger.hpp"

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
typedef ::World<WorldTraits> World;
typedef ::Model Model;
typedef ::ParticleModel ParticleModel;
typedef ::NetworkRules NetworkRules;
typedef NetworkRules::reaction_rule_generator ReactionRuleGenerator;
typedef World::transaction_type Transaction;
typedef World::base_type::base_type ParticleContainer; 
typedef ::TransactionImpl<ParticleContainer> TransactionImpl;
typedef ::EGFRDSimulatorTraitsBase<World> EGFRDSimulatorTraits;
typedef ::ParticleSimulator<EGFRDSimulatorTraits> ParticleSimulator;
typedef ::BDSimulator<EGFRDSimulatorTraits> BDSimulator;
typedef ::EGFRDSimulator<EGFRDSimulatorTraits> EGFRDSimulator;
typedef ::MultiParticleContainer<EGFRDSimulatorTraits> MultiParticleContainer;
typedef EGFRDSimulator::box_type Box;
typedef EGFRDSimulator::sphere_type Sphere;
typedef EGFRDSimulator::cylinder_type Cylinder;
typedef EGFRDSimulator::plane_type Plane;
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
typedef ::StructureUtils<EGFRDSimulator> StructureUtils;
typedef EGFRDSimulator::particle_simulation_structure_type ParticleSimulationStructure;
typedef EGFRDSimulator::surface_type Surface;
typedef EGFRDSimulator::region_type Region;
typedef EGFRDSimulator::planar_surface_type PlanarSurface;
typedef EGFRDSimulator::spherical_surface_type SphericalSurface;
typedef EGFRDSimulator::cylindrical_surface_type CylindricalSurface;
typedef EGFRDSimulator::cuboidal_region_type CuboidalRegion;
typedef EGFRDSimulatorTraits::reaction_record_type ReactionRecord;
typedef EGFRDSimulatorTraits::reaction_recorder_type ReactionRecorder;
typedef EGFRDSimulatorTraits::volume_clearer_type VolumeClearer;
typedef ::Logger Logger;
typedef ::LogAppender LogAppender;
typedef ::LoggerManager LoggerManager;

} // namespace binding

#endif /* TRAITS_HPP */
