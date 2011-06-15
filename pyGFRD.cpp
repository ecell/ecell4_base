#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <exception>
#include <stdexcept>

//#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <boost/lexical_cast.hpp>
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>
#include <boost/python/refcount.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/reference_existing_object.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/converter/object_manager.hpp>

#include <numpy/arrayobject.h>

#include "utils.hpp"
#include "peer/utils.hpp"
#include "peer/converters/tuple.hpp"
#include "peer/numpy/wrapped_multi_array.hpp"
#include "peer/numpy/scalar_converters.hpp"

#include "binding/bd_propagator_class.hpp"
#include "binding/binding_common.hpp"
#include "binding/box_class.hpp"
#include "binding/cylinder_class.hpp"
#include "binding/domain_id_class.hpp"
#include "binding/domain_classes.hpp"
#include "binding/egfrd_simulator_classes.hpp"
#include "binding/bd_simulator_classes.hpp"
#include "binding/event_classes.hpp"
#include "binding/exception_classes.hpp"
#include "binding/matrix_space_classes.hpp"
#include "binding/model_class.hpp"
#include "binding/module_functions.hpp"
#include "binding/multi_particle_container_class.hpp"
#include "binding/network_rules_class.hpp"
#include "binding/network_rules_wrapper_class.hpp"
#include "binding/particle_class.hpp"
#include "binding/particle_container_class.hpp"
#include "binding/particle_id_class.hpp"
#include "binding/particle_model_class.hpp"
#include "binding/particle_simulator_classes.hpp"
#include "binding/plane_class.hpp"
#include "binding/position_converters.hpp"
#include "binding/python_logger_classes.hpp"
#include "binding/random_number_generator_class.hpp"
#include "binding/volume_clearer_classes.hpp"
#include "binding/reaction_record_classes.hpp"
#include "binding/shape_converters.hpp"
#include "binding/shell_classes.hpp"
#include "binding/shell_id_class.hpp"
#include "binding/species_id_class.hpp"
#include "binding/species_type_class.hpp"
#include "binding/sphere_class.hpp"
#include "binding/structure_classes.hpp"
#include "binding/transaction_classes.hpp"
#include "binding/world_class.hpp"

namespace b = binding;

BOOST_PYTHON_MODULE(_gfrd)
{
    using namespace boost::python;

    import_array();

    // GSL error handler: is this the best place for this?
    gsl_set_error_handler(&gsl_error_handler);

    peer::util::register_std_exception_translator();
    peer::util::register_exception_translator<PyExc_IndexError, 
        std::out_of_range>();

    b::register_model_class();
    b::register_bd_propagator_class();
    b::register_box_class();
    b::register_domain_id_class();
    b::register_domain_classes();
    b::register_exception_classes();
    b::register_spherical_shell_container_class();
    b::register_plane_class();
    b::register_cylinder_class();
    b::register_cylindrical_shell_container_class();
    b::register_network_rules_class();
    b::register_network_rules_wrapper_class();
    b::register_particle_class();
    b::register_particle_id_class();
    b::register_particle_model_class();
    b::register_position_converters();
    b::register_event_class();
    b::register_python_event_class();
    b::register_event_scheduler_class();
    b::register_random_number_generator_class();
    b::register_sphere_class();
    b::register_sphere_converters();
    b::register_spherical_shell_class();
    b::register_cylindrical_shell_class();
    b::register_shell_id_class();
    b::register_species_id_class();
    b::register_species_type_class();
    b::register_particle_container_class();
    b::register_multi_particle_container_class();
    b::register_transaction_classes();
    b::register_world_class();
    b::register_structure_classes();
    b::register_module_functions();
    b::register_volume_clearer_classes();
    b::register_reaction_record_classes();
    b::register_particle_simulator_classes();
    b::register_egfrd_simulator_classes();
    b::register_bd_simulator_classes();
    b::register_python_logger_classes();

    peer::util::register_seq_wrapped_multi_array_converter<b::Length>();
    peer::util::register_ndarray_wrapped_multi_array_converter<b::Length, 2>();
    peer::util::register_ndarray_wrapped_multi_array_converter<b::Length, 3>();

    peer::util::register_scalar_to_native_converters();
}
