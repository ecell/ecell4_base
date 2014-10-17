#ifndef __ECELL4_EGFRD_EGFRD_HPP
#define __ECELL4_EGFRD_EGFRD_HPP

#include <ecell4/core/types.hpp>

#include "World.hpp"
#include "EGFRDSimulator.hpp"

namespace ecell4
{

namespace egfrd
{

typedef ::World< ::CyclicWorldTraits<Real> > EGFRDWorld;
typedef EGFRDWorld::molecule_info_type MoleculeInfo;
typedef ::EGFRDSimulator< ::EGFRDSimulatorTraitsBase<EGFRDWorld> > EGFRDSimulator;

} // egfrd

} // ecell4

#endif /* __ECELL4_EGFRD_EGFRD_HPP */
