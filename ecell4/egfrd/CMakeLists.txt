set(CPP_FILES
    Logger.cpp ConsoleAppender.cpp utils.cpp findRoot.cpp funcSum.cpp
    GreensFunction1DAbsAbs.cpp GreensFunction1DRadAbs.cpp
    PairGreensFunction.cpp GreensFunction3D.cpp GreensFunction3DAbs.cpp
    GreensFunction3DAbsSym.cpp GreensFunction3DRadAbsBase.cpp
    GreensFunction3DRadAbs.cpp GreensFunction3DRadInf.cpp
    GreensFunction3DSym.cpp SphericalBesselGenerator.cpp CylindricalBesselGenerator.cpp)

set(HPP_FILES
    Multi.hpp ReactionRecord.hpp sorted_list.hpp factorial.hpp ParticleContainerBase.hpp NetworkRulesAdapter.hpp Logger.hpp MatrixSpace.hpp DomainFactory.hpp abstract_set.hpp twofold_container.hpp ParticleTraits.hpp DomainID.hpp geometry.hpp generator.hpp ReactionRecorder.hpp AnalyticalPair.hpp utils.hpp DomainUtils.hpp EGFRDSimulator.hpp BDSimulator.hpp freeFunctions.hpp Real3Type.hpp ReactionRuleInfo.hpp BDPropagator.hpp VolumeClearer.hpp funcSum.hpp World.hpp OldDefs.hpp AnalyticalSingle.hpp ShapedDomain.hpp ParticleContainer.hpp ParticleSimulator.hpp Defs.hpp Pair.hpp Domain.hpp filters.hpp ShellID.hpp Shell.hpp findRoot.hpp ConsoleAppender.hpp Single.hpp Transaction.hpp linear_algebra.hpp exceptions.hpp ReactionRecorderWrapper.hpp
    GreensFunction3D.hpp GreensFunction3DSym.hpp GreensFunction3DRadInf.hpp SphericalBesselGenerator.hpp GreensFunction1DAbsAbs.hpp PairGreensFunction.hpp GreensFunction3DRadAbsBase.hpp GreensFunction3DAbs.hpp GreensFunction3DRadAbs.hpp GreensFunction3DAbsSym.hpp GreensFunction1DRadAbs.hpp GreensFunction.hpp CylindricalBesselGenerator.hpp HalfOrderBesselGenerator.hpp
    compat.h
    egfrd.hpp structures.hpp)

set(UTILS_HPP_FILES
    utils/base_type_walker.hpp utils/memberwise_compare.hpp utils/array_helper.hpp utils/range_support.hpp utils/random.hpp utils/pair.hpp utils/fun_wrappers.hpp utils/reset.hpp utils/math.hpp utils/swap.hpp utils/unassignable_adapter.hpp utils/stringizer.hpp utils/array_traits.hpp utils/fun_composition.hpp utils/range.hpp utils/assoc_container_traits.hpp utils/reference_or_instance.hpp utils/pointer_as_ref.hpp utils/get_default_impl.hpp utils/map_adapter.hpp utils/pointer_preds.hpp)

# CPP_FILES BasicNetworkRulesImpl.cpp NetworkRules.cpp SpeciesType.cpp StructureType.cpp
# HPP_FILES StructureUtils.hpp EGFRDSimulatorFactory.hpp ParticleSimulatorFactory.hpp Model.hpp ParticleModel.hpp BasicNetworkRulesImpl.hpp NetworkRules.hpp ReactionRule.hpp SpeciesInfo.hpp BDSimulator.hpp SpeciesTypeID.hpp SpeciesType.hpp ParticleID.hpp StructureType.hpp bessel.hpp PyEventScheduler.hpp NetworkRulesWrapper.hpp BasicReactionRecorderImpl.hpp DynamicPriorityQueue.hpp EventScheduler.hpp SerialIDGenerator.hpp Surface.hpp SphericalSurface.hpp Structure.hpp CylindricalSurface.hpp PlanarSurface.hpp CuboidalRegion.hpp ParticleSimulationStructure.hpp Region.hpp Box.hpp Plane.hpp Cylinder.hpp Shape.hpp Sphere.hpp Point.hpp Vector3.hpp Identifier.hpp

# find_library(BOOST_REGEX_LIBRARIES boost_regex)

check_cxx_source_compiles("#include <math.h>
int main() { double a = INFINITY; return (int)a * 0; }"
    HAVE_DECL_INFINITY)
check_cxx_source_compiles("#include <math.h>
int main() { isfinite(0.); return 0; }"
    HAVE_ISFINITE)
check_cxx_source_compiles("#include <math.h>
int main() { double a, b; sincos(0., &a, &b); return 0; }"
    HAVE_SINCOS)

configure_file(
    ${PROJECT_SOURCE_DIR}/ecell4/egfrd/config.h.in
    ${PROJECT_BINARY_DIR}/ecell4/egfrd/config.h)

execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/ecell4/egfrd)

if (NOT NO_BESSEL_TABLE)
add_subdirectory(tablegen)
endif()

include_directories(${PROJECT_SOURCE_DIR}/ecell4/egfrd ${PROJECT_BINARY_DIR}/ecell4/egfrd)

# add_definitions(-DHAVE_CONFIG_H -DHAVE_INLINE)
add_definitions(-DHAVE_CONFIG_H)

add_library(ecell4-egfrd SHARED
    ${CPP_FILES} ${HPP_FILES} ${UTILS_HPP_FILES}
    "${PROJECT_BINARY_DIR}/ecell4/egfrd/config.h")
target_link_libraries(ecell4-egfrd ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} ecell4-core)
# target_link_libraries(ecell4-egfrd ecell4-core ${BOOST_REGEX_LIBRARIES})

add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-egfrd DESTINATION lib)

install(FILES ${HPP_FILES} "${PROJECT_BINARY_DIR}/ecell4/egfrd/config.h"
    DESTINATION "include/ecell4/egfrd")

if (NOT NO_BESSEL_TABLE)
add_dependencies(ecell4-egfrd BesselTables)

install(FILES ${HPP_FILES}
    "${PROJECT_BINARY_DIR}/ecell4/egfrd/SphericalBesselTable.hpp"
    "${PROJECT_BINARY_DIR}/ecell4/egfrd/CylindricalBesselTable.hpp"
    DESTINATION "include/ecell4/egfrd")
endif()

install(FILES ${UTILS_HPP_FILES} DESTINATION "include/ecell4/egfrd/utils")
