set(CPP_FILES
    Logger.cpp ConsoleAppender.cpp utils.cpp findRoot.cpp
    funcSum.cpp GreensFunction1DAbsAbs.cpp GreensFunction1DRadAbs.cpp
    PairGreensFunction.cpp GreensFunction3D.cpp GreensFunction3DAbs.cpp
    GreensFunction3DAbsSym.cpp GreensFunction3DRadAbsBase.cpp
    GreensFunction3DRadAbs.cpp GreensFunction3DRadInf.cpp
    GreensFunction3DSym.cpp SphericalBesselGenerator.cpp
    BasicNetworkRulesImpl.cpp NetworkRules.cpp SpeciesType.cpp
    StructureType.cpp)

set(HPP_FILES
    Surface.hpp Multi.hpp Box.hpp Plane.hpp Vector3.hpp ReactionRecord.hpp ParticleSimulatorFactory.hpp sorted_list.hpp Model.hpp factorial.hpp ParticleContainerBase.hpp NetworkRulesAdapter.hpp Logger.hpp MatrixSpace.hpp SphericalSurface.hpp DomainFactory.hpp abstract_set.hpp BasicReactionRecorderImpl.hpp SpeciesTypeID.hpp twofold_container.hpp ParticleTraits.hpp DomainID.hpp GreensFunction3D.hpp BDSimulator.hpp geometry.hpp generator.hpp Cylinder.hpp ReactionRecorder.hpp GreensFunction1DAbsAbs.hpp PairGreensFunction.hpp AnalyticalPair.hpp utils.hpp SerialIDGenerator.hpp DomainUtils.hpp EGFRDSimulator.hpp freeFunctions.hpp GreensFunction3DRadAbsBase.hpp Position3Type.hpp Structure.hpp ReactionRuleInfo.hpp CylindricalSurface.hpp BDPropagator.hpp VolumeClearer.hpp funcSum.hpp PyEventScheduler.hpp StructureUtils.hpp SpeciesInfo.hpp GreensFunction3DAbs.hpp CylindricalBesselGenerator.hpp World.hpp GreensFunction3DRadAbs.hpp Identifier.hpp GreensFunction3DAbsSym.hpp Point.hpp PlanarSurface.hpp GreensFunction1DRadAbs.hpp OldDefs.hpp CuboidalRegion.hpp AnalyticalSingle.hpp HalfOrderBesselGenerator.hpp EventScheduler.hpp bessel.hpp EGFRDSimulatorFactory.hpp BasicNetworkRulesImpl.hpp ShapedDomain.hpp ParticleContainer.hpp ParticleSimulationStructure.hpp ParticleModel.hpp ParticleSimulator.hpp Region.hpp ParticleID.hpp SpeciesType.hpp Defs.hpp Pair.hpp Domain.hpp GreensFunction3DSym.hpp NetworkRulesWrapper.hpp DynamicPriorityQueue.hpp filters.hpp ShellID.hpp ReactionRule.hpp Shape.hpp Shell.hpp findRoot.hpp ConsoleAppender.hpp GreensFunction.hpp Single.hpp GreensFunction3DRadInf.hpp StructureType.hpp NetworkRules.hpp Transaction.hpp linear_algebra.hpp Sphere.hpp exceptions.hpp SphericalBesselGenerator.hpp
    utils/base_type_walker.hpp utils/memberwise_compare.hpp utils/array_helper.hpp utils/range_support.hpp utils/random.hpp utils/pair.hpp utils/fun_wrappers.hpp utils/reset.hpp utils/math.hpp utils/swap.hpp utils/unassignable_adapter.hpp utils/stringizer.hpp utils/array_traits.hpp utils/fun_composition.hpp utils/range.hpp utils/assoc_container_traits.hpp utils/reference_or_instance.hpp utils/pointer_as_ref.hpp utils/get_default_impl.hpp utils/map_adapter.hpp utils/pointer_preds.hpp
    compat.h)

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
    ${CPP_FILES} ${HPP_FILES} "${PROJECT_BINARY_DIR}/ecell4/epdp/config.h")
target_link_libraries(ecell4-epdp ecell4-core ${BOOST_REGEX_LIBRARIES})

# add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-epdp DESTINATION lib)
add_dependencies(ecell4-epdp BesselTables)

install(FILES ${HPP_FILES} "${PROJECT_BINARY_DIR}/ecell4/epdp/config.h"
    "${PROJECT_BINARY_DIR}/ecell4/epdp/SphericalBesselTable.hpp"
    "${PROJECT_BINARY_DIR}/ecell4/epdp/CylindricalBesselTable.hpp"
    DESTINATION "include/ecell4/epdp")
