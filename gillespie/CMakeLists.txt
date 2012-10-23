CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(gillespie_solver_unit)

SET(BUILD_GILLESPIE_UNIT ON)
SET(CMAKE_CXX_FLAGS_DEBUG "-g -pg")

# build libgillespie_solver_unit 
SET(serial "1.0.0")
SET(soserial "1")
ADD_LIBRARY(gillespie_solver SHARED  Gillespie.cpp)
SET_TARGET_PROPERTIES(gillespie_solver PROPERTIES VERSION ${serial} SOVERSION ${soserial})
TARGET_LINK_LIBRARIES(gillespie_solver gsl gslcblas)
INSTALL(TARGETS gillespie_solver LIBRARY DESTINATION lib)
INSTALL(FILES Gillespie.hpp DESTINATION include)

# build small binary using libgillespie.so or dylib(mac)
if (BUILD_GILLESPIE_UNIT)
	ADD_EXECUTABLE(gillespie_unit_test ./test.cpp)
	TARGET_LINK_LIBRARIES(gillespie_unit_test gillespie_solver)
endif()
