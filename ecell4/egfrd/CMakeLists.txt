execute_process(COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/ecell4/egfrd)

if (NO_SHARED)
    return()
endif()

set(CPP_FILES
    Logger.cpp ConsoleAppender.cpp utils.cpp
    )

set(HPP_FILES
    Multi.hpp ReactionRecord.hpp sorted_list.hpp factorial.hpp ParticleContainerBase.hpp NetworkRulesAdapter.hpp Logger.hpp MatrixSpace.hpp DomainFactory.hpp abstract_set.hpp twofold_container.hpp ParticleTraits.hpp DomainID.hpp geometry.hpp generator.hpp ReactionRecorder.hpp AnalyticalPair.hpp utils.hpp DomainUtils.hpp EGFRDSimulator.hpp BDSimulator.hpp Real3Type.hpp ReactionRuleInfo.hpp BDPropagator.hpp VolumeClearer.hpp World.hpp AnalyticalSingle.hpp ShapedDomain.hpp ParticleContainer.hpp ParticleSimulator.hpp Defs.hpp Pair.hpp Domain.hpp filters.hpp ShellID.hpp Shell.hpp ConsoleAppender.hpp Single.hpp Transaction.hpp linear_algebra.hpp exceptions.hpp ReactionRecorderWrapper.hpp
    egfrd.hpp structures.hpp)

set(UTILS_HPP_FILES
    utils/base_type_walker.hpp utils/memberwise_compare.hpp utils/array_helper.hpp utils/range_support.hpp utils/random.hpp utils/pair.hpp utils/fun_wrappers.hpp utils/reset.hpp utils/math.hpp utils/swap.hpp utils/unassignable_adapter.hpp utils/stringizer.hpp utils/array_traits.hpp utils/fun_composition.hpp utils/range.hpp utils/assoc_container_traits.hpp utils/reference_or_instance.hpp utils/pointer_as_ref.hpp utils/get_default_impl.hpp utils/map_adapter.hpp utils/pointer_preds.hpp)

include_directories(${PROJECT_SOURCE_DIR}/ecell4/egfrd ${PROJECT_BINARY_DIR}/ecell4/egfrd)

add_library(ecell4-egfrd SHARED
    ${CPP_FILES} ${HPP_FILES} ${UTILS_HPP_FILES})
target_link_libraries(ecell4-egfrd ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES} ecell4-core greens_functions)

set(ECELL4_SHARED_DIRS ${CMAKE_CURRENT_BINARY_DIR}:${ECELL4_SHARED_DIRS} PARENT_SCOPE)

# add_subdirectory(tests)
add_subdirectory(samples)

install(TARGETS ecell4-egfrd DESTINATION lib)

install(FILES ${HPP_FILES} DESTINATION "include/ecell4/egfrd")

install(FILES ${UTILS_HPP_FILES} DESTINATION "include/ecell4/egfrd/utils")
