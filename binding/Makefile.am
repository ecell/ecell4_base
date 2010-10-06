NUMPY_INCLUDE_DIR = @NUMPY_INCLUDE_DIR@
PYTHON_INCLUDES = @PYTHON_INCLUDES@

INCLUDES = ${PYTHON_INCLUDES} -I${NUMPY_INCLUDE_DIR}
AM_CXXFLAGS = @CXXFLAGS@ @GSL_CFLAGS@ -I.. -Wall -g -Wstrict-aliasing=0 -Wno-invalid-offsetof
AM_CPPFLAGS = @CPPFLAGS@ -DNO_IMPORT

noinst_LTLIBRARIES = libbinding_utils.la

noinst_HEADERS = \
	BDPropagator.hpp \
	bd_propagator_class.hpp \
	Box.hpp \
	box_class.hpp \
	binding_common.hpp \
	Cylinder.hpp \
	cylinder_class.hpp \
	CylindricalSurface.hpp \
	CuboidalRegion.hpp \
	domain_classes.hpp \
	domain_id_class.hpp \
	exception_classes.hpp \
	Identifier.hpp \
	matrix_space_classes.hpp \
	MatrixSpace.hpp \
	Model.hpp \
	model_class.hpp \
	module_functions.hpp \
	MultiParticleContainer.hpp \
	multi_particle_container_class.hpp \
	network_rules_class.hpp \
	NetworkRules.hpp \
	network_rules_wrapper_class.hpp \
	NetworkRulesWrapper.hpp \
	Pair.hpp \
	particle_class.hpp \
	ParticleContainer.hpp \
	Particle.hpp \
	particle_id_class.hpp \
	particle_container_class.hpp \
	particle_model_class.hpp \
	position_converters.hpp \
	PlanarSurface.hpp \
	plane_class.hpp \
	event_classes.hpp \
	Event.hpp \
	EventScheduler.hpp \
	random_number_generator_class.hpp \
	RandomNumberGenerator.hpp \
	ReactionRecord.hpp \
	reaction_recorder_converter.hpp \
	reaction_record_classes.hpp \
	ReactionRule.hpp \
	ReactionRuleInfo.hpp \
	SerialIDGenerator.hpp \
	ShapedDomain.hpp
	shape_converters.hpp \
	shell_classes.hpp \
	Shell.hpp \
	shell_id_class.hpp \
	Single.hpp \
	species_id_class.hpp \
	SpeciesInfo.hpp \
	species_type_class.hpp \
	SpeciesType.hpp \
	Sphere.hpp \
	sphere_class.hpp \
	SphericalSurface.hpp \
	StructureType.hpp \
	structure_classes.hpp \
	structure_type_class.hpp \
	transaction_classes.hpp \
	Transaction.hpp \
	world_class.hpp \
	World.hpp

libbinding_utils_la_SOURCES = \
	bd_propagator_class.cpp \
	box_class.cpp \
	cylinder_class.cpp \
	domain_classes.cpp \
	domain_id_class.cpp \
	egfrd_simulator_classes.cpp \
	exception_classes.cpp \
	matrix_space_classes.cpp \
	model_class.cpp \
	module_functions.cpp \
	multi_particle_container_class.cpp \
	network_rules_class.cpp \
	network_rules_wrapper_class.cpp \
	particle_class.cpp \
	particle_container_class.cpp \
	particle_id_class.cpp \
	particle_model_class.cpp \
	plane_class.cpp \
	position_converters.cpp \
	event_classes.cpp \
	random_number_generator_class.cpp \
	reaction_record_classes.cpp \
	shape_converters.cpp \
	shell_classes.cpp \
	shell_id_class.cpp \
	species_id_class.cpp \
	species_type_class.cpp \
	sphere_class.cpp \
	structure_classes.cpp \
	structure_type_class.cpp \
	transaction_classes.cpp \
	world_class.cpp
