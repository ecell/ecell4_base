set(CPP_FILES
    Logger.cpp ConsoleAppender.cpp utils.cpp findRoot.cpp funcSum.cpp
    GreensFunction1DAbsAbs.cpp GreensFunction1DRadAbs.cpp
    PairGreensFunction.cpp GreensFunction3D.cpp GreensFunction3DAbs.cpp
    GreensFunction3DAbsSym.cpp GreensFunction3DRadAbsBase.cpp
    GreensFunction3DRadAbs.cpp GreensFunction3DRadInf.cpp
    GreensFunction3DSym.cpp SphericalBesselGenerator.cpp)

set(HPP_FILES
    Surface.hpp Multi.hpp Box.hpp Plane.hpp Vector3.hpp ReactionRecord.hpp sorted_list.hpp factorial.hpp ParticleContainerBase.hpp NetworkRulesAdapter.hpp Logger.hpp MatrixSpace.hpp SphericalSurface.hpp DomainFactory.hpp abstract_set.hpp twofold_container.hpp ParticleTraits.hpp DomainID.hpp geometry.hpp generator.hpp Cylinder.hpp ReactionRecorder.hpp AnalyticalPair.hpp utils.hpp DomainUtils.hpp EGFRDSimulator.hpp freeFunctions.hpp Position3Type.hpp Structure.hpp ReactionRuleInfo.hpp CylindricalSurface.hpp BDPropagator.hpp VolumeClearer.hpp funcSum.hpp World.hpp Identifier.hpp Point.hpp PlanarSurface.hpp OldDefs.hpp CuboidalRegion.hpp AnalyticalSingle.hpp ShapedDomain.hpp ParticleContainer.hpp ParticleSimulationStructure.hpp ParticleSimulator.hpp Region.hpp Defs.hpp Pair.hpp Domain.hpp filters.hpp ShellID.hpp Shape.hpp Shell.hpp findRoot.hpp ConsoleAppender.hpp Single.hpp Transaction.hpp linear_algebra.hpp Sphere.hpp exceptions.hpp ReactionRecorderWrapper.hpp
    GreensFunction3D.hpp GreensFunction3DSym.hpp GreensFunction3DRadInf.hpp SphericalBesselGenerator.hpp GreensFunction1DAbsAbs.hpp PairGreensFunction.hpp GreensFunction3DRadAbsBase.hpp GreensFunction3DAbs.hpp GreensFunction3DRadAbs.hpp GreensFunction3DAbsSym.hpp GreensFunction1DRadAbs.hpp GreensFunction.hpp CylindricalBesselGenerator.hpp HalfOrderBesselGenerator.hpp
    compat.h
    egfrd.hpp)

set(UTILS_HPP_FILES
    utils/base_type_walker.hpp utils/memberwise_compare.hpp utils/array_helper.hpp utils/range_support.hpp utils/random.hpp utils/pair.hpp utils/fun_wrappers.hpp utils/reset.hpp utils/math.hpp utils/swap.hpp utils/unassignable_adapter.hpp utils/stringizer.hpp utils/array_traits.hpp utils/fun_composition.hpp utils/range.hpp utils/assoc_container_traits.hpp utils/reference_or_instance.hpp utils/pointer_as_ref.hpp utils/get_default_impl.hpp utils/map_adapter.hpp utils/pointer_preds.hpp)

# CPP_FILES BasicNetworkRulesImpl.cpp NetworkRules.cpp SpeciesType.cpp StructureType.cpp
# HPP_FILES StructureUtils.hpp EGFRDSimulatorFactory.hpp ParticleSimulatorFactory.hpp Model.hpp ParticleModel.hpp BasicNetworkRulesImpl.hpp NetworkRules.hpp ReactionRule.hpp SpeciesInfo.hpp BDSimulator.hpp SpeciesTypeID.hpp SpeciesType.hpp ParticleID.hpp StructureType.hpp bessel.hpp PyEventScheduler.hpp NetworkRulesWrapper.hpp BasicReactionRecorderImpl.hpp DynamicPriorityQueue.hpp EventScheduler.hpp SerialIDGenerator.hpp

find_library(BOOST_REGEX_LIBRARIES boost_regex)

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
    ${PROJECT_SOURCE_DIR}/ecell4/epdp/config.h.in
    ${PROJECT_BINARY_DIR}/ecell4/epdp/config.h)

execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/ecell4/epdp/epdp)

add_custom_command(
    OUTPUT SphericalBesselTable.hpp
    COMMAND python make_sjy_table.py ${PROJECT_BINARY_DIR}/ecell4/epdp/SphericalBesselTable.hpp)
add_custom_command(
    OUTPUT CylindricalBesselTable.hpp
    COMMAND python make_cjy_table.py ${PROJECT_BINARY_DIR}/ecell4/epdp/CylindricalBesselTable.hpp)
add_custom_target(BesselTables ALL DEPENDS ${PROJECT_BINARY_DIR}/ecell4/epdp/SphericalBesselTable.hpp ${PROJECT_BINARY_DIR}/ecell4/epdp/CylindricalBesselTable.hpp)

include_directories(${PROJECT_SOURCE_DIR}/ecell4/epdp ${PROJECT_BINARY_DIR}/ecell4/epdp)

add_definitions(-DHAVE_CONFIG_H -DHAVE_INLINE)

add_library(ecell4-epdp SHARED
    ${CPP_FILES} ${HPP_FILES} ${UTILS_HPP_FILES}
    "${PROJECT_BINARY_DIR}/ecell4/epdp/config.h")
target_link_libraries(ecell4-epdp ecell4-core ${BOOST_REGEX_LIBRARIES})

# add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-epdp DESTINATION lib)
add_dependencies(ecell4-epdp BesselTables)

install(FILES ${HPP_FILES} "${PROJECT_BINARY_DIR}/ecell4/epdp/config.h"
    "${PROJECT_BINARY_DIR}/ecell4/epdp/SphericalBesselTable.hpp"
    "${PROJECT_BINARY_DIR}/ecell4/epdp/CylindricalBesselTable.hpp"
    DESTINATION "include/ecell4/epdp")
install(FILES ${UTILS_HPP_FILES} DESTINATION "include/ecell4/epdp/utils")
