#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <exception>
#include <stdexcept>

//#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <boost/lexical_cast.hpp>
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>
#include <boost/python/refcount.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>

#include <numpy/arrayobject.h>

#include "PyEventScheduler.hpp"
#include "freeFunctions.hpp"
#include "FreeGreensFunction.hpp"
#include "FirstPassageGreensFunction.hpp"
#include "BasicPairGreensFunction.hpp"
#include "FreePairGreensFunction.hpp"
#include "FirstPassagePairGreensFunction.hpp"
#include "FirstPassageNoCollisionPairGreensFunction.hpp"

#include "array_traits.hpp"
#include "Vector3.hpp"
#include "utils.hpp"
#include "Model.hpp"

#include "peer/utils.hpp"
#include "peer/tuple_converters.hpp"
#include "peer/numpy/wrapped_multi_array.hpp"
#include "peer/MatrixSpace.hpp"
#include "peer/SpeciesType.hpp"
#include "peer/Identifier.hpp"
#include "peer/ReactionRule.hpp"
#include "peer/GeneratorIteratorWrapper.hpp"

typedef Real length_type;
typedef Vector3<length_type> position_type;

struct position_to_ndarray_converter
{
    typedef position_type native_type;
    
    static PyObject* convert( const native_type& p )
    {
        static const npy_intp dims[1] = { native_type::size() };
        void* data( PyDataMem_NEW( native_type::size() * sizeof( position_type::value_type ) ) );
        memcpy( data, static_cast<const void*>( p.data() ),
                native_type::size() * sizeof( position_type::value_type ) );
        PyObject* array( PyArray_New( &PyArray_Type, 1, 
                                      const_cast<npy_intp*>( dims ),
                                      peer::util::get_numpy_typecode<
											position_type::value_type >
                                      ::value, NULL,
                                      data, 0, NPY_CARRAY, NULL ) );
        reinterpret_cast<PyArrayObject*>( array )->flags |= NPY_OWNDATA;
        return array;
    }
};

struct ndarray_to_position_converter
{
    typedef position_type native_type;
    
    static void* convertible(PyObject* ptr)
    {
        if (!PyArray_Check(ptr))
        {
            return NULL;
        }
        
        PyObject* retval(PyArray_CastToType(
                             reinterpret_cast<PyArrayObject*>(ptr),
                             PyArray_DescrFromType(
                                 peer::util::get_numpy_typecode<
										position_type::value_type >::value ),
                                 0) );
        if (!retval)
        {
            return NULL;
        }
        
        if (PyObject_Size(retval) != 3)
        {
            boost::python::decref(retval);
            return NULL;
        }
        
        return retval;
    }
    
    static void construct(PyObject* ptr,
                          boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        PyArrayObject* array_obj = static_cast<PyArrayObject*>(
            data->stage1.convertible);
        data->stage1.convertible = new(data->storage.bytes) native_type(
            reinterpret_cast<position_type::value_type*>(PyArray_DATA(array_obj)));
        boost::python::decref(reinterpret_cast<PyObject*>(array_obj));
    }
};

struct seq_to_position_converter
{
    typedef position_type native_type;
    
    static void* convertible(PyObject* ptr)
    {
        if (!PySequence_Check(ptr))
        {
            return NULL;
        }
        
        if (PySequence_Size(ptr) != 3)
        {
            return NULL;
        }
        
        return ptr;
    }
    
    static void construct(PyObject* ptr,
                          boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        data->stage1.convertible = new(data->storage.bytes) native_type(
            PyFloat_AsDouble(PySequence_GetItem(ptr, 0)),
            PyFloat_AsDouble(PySequence_GetItem(ptr, 1)),
            PyFloat_AsDouble(PySequence_GetItem(ptr, 2)));
    }
};

struct sphere_to_python_converter
{
    typedef Sphere<length_type> native_type;

    static PyObject* convert(native_type const& v)
    {
        return boost::python::incref(
            boost::python::object(boost::make_tuple(
                v.position(), v.radius())).ptr());
    }
};

struct python_to_sphere_converter
{
    typedef Sphere<length_type> native_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PyTuple_Check(pyo))
        {
            return 0;
        }
        if (PyTuple_Size(pyo) != 2)
        {
            return 0;
        }
        return pyo;
    }

    static void construct(PyObject* pyo, 
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        PyObject* items[] = { PyTuple_GetItem(pyo, 0), PyTuple_GetItem(pyo, 1) };
        void* storage(reinterpret_cast<
            boost::python::converter::rvalue_from_python_storage<
                native_type >*
            >(data)->storage.bytes);
        new (storage) native_type (
            boost::python::extract<
                native_type::position_type>(items[0]),
            PyFloat_AsDouble(items[1]));
        data->convertible = storage;
    }
};


BOOST_PYTHON_MODULE( _gfrd )
{
    using namespace boost::python;
    typedef Real length_type;
    typedef Vector3< length_type > vector_type;

    import_array();

    // GSL error handler: is this the best place for this?
    gsl_set_error_handler( &gsl_error_handler );

  
    peer::util::register_std_exception_translator();

    peer::util::register_seq_wrapped_multi_array_converter<length_type>();
    // peer::util::register_ndarray_wrapped_multi_array_converter<length_type, 1>();
    peer::util::register_ndarray_wrapped_multi_array_converter<length_type, 2>();
    peer::util::register_ndarray_wrapped_multi_array_converter<length_type, 3>();

    peer::util::register_tuple_converter< boost::tuple< Real, EventType > >();

    // free functions
    def( "p_irr", p_irr );
    def( "p_survival_irr", p_survival_irr );
    def( "p_theta_free", p_theta_free );
    def( "ip_theta_free", ip_theta_free );
    def( "g_bd", g_bd );
    def( "I_bd", I_bd );
    def( "I_bd_r", I_bd_r );
    def( "drawR_gbd", drawR_gbd );
    def( "p_reaction_irr", __p_reaction_irr );
    def( "p_reaction_irr_t_inf", __p_reaction_irr_t_inf );

    class_<PyEvent>( "PyEvent", 
                     init<const Real, const object&, const object&>() )
        .def( "setTime", &PyEvent::setTime )
        .def( "getTime", &PyEvent::getTime )
        .def( "getObj", &PyEvent::getObj,
              return_value_policy<copy_const_reference>() )
        .def( "getArg", &PyEvent::getArg,
              return_value_policy<copy_const_reference>() )
        ;


    typedef const PyEventScheduler::Event& 
        (PyEventScheduler::*geteventrefsig)() const;
    typedef const PyEventScheduler::Event& 
        (PyEventScheduler::*geteventrefbyindexsig)
        ( const PyEventScheduler::EventID );

    class_<PyEventScheduler, boost::noncopyable>( "EventScheduler" )
        .def( "getTime", &PyEventScheduler::getTime )
        .def( "getTopTime", &PyEventScheduler::getTopTime )
        .def( "getSize", &PyEventScheduler::getSize )
        .def( "getTopEvent", geteventrefsig( &PyEventScheduler::getTopEvent ),
              return_value_policy<copy_const_reference>() )
        .def( "getTopID", &PyEventScheduler::getTopID )
        .def( "peekSecondEvent", 
              geteventrefsig( &PyEventScheduler::peekSecondEvent ),
              return_value_policy<copy_const_reference>() )
        .def( "getEvent", geteventrefbyindexsig( &PyEventScheduler::getEvent ),
              return_value_policy<copy_const_reference>() )
        .def( "getEventByIndex", &PyEventScheduler::getEventByIndex,
              return_value_policy<copy_const_reference>() )
        .def( "step", &PyEventScheduler::step )
        .def( "clear", &PyEventScheduler::clear )
        .def( "addEvent", &PyEventScheduler::addEvent )
        .def( "removeEvent", &PyEventScheduler::removeEvent )
        .def( "updateEventTime", &PyEventScheduler::updateEventTime )
//        .def( "updateAllEventDependency", 
//              &PyEventScheduler::updateAllEventDependency )
        .def( "check", &PyEventScheduler::check )
        ;


    class_<FreeGreensFunction>( "FreeGreensFunction",
                                init<const Real>() )
        .def( "getD", &FreeGreensFunction::getD )
        .def( "seta", &FreeGreensFunction::seta )
        .def( "drawTime", &FreeGreensFunction::drawTime )
        .def( "drawR", &FreeGreensFunction::drawR )
        .def( "p_r", &FreeGreensFunction::p_r )
        .def( "ip_r", &FreeGreensFunction::ip_r )
        .def( "dump", &FreeGreensFunction::dump )
        ;

    class_<FirstPassageGreensFunction>( "FirstPassageGreensFunction",
                                        init<const Real>() )
        .def( "getD", &FirstPassageGreensFunction::getD )
        .def( "seta", &FirstPassageGreensFunction::seta )
        .def( "geta", &FirstPassageGreensFunction::geta )
        .def( "drawTime", &FirstPassageGreensFunction::drawTime )
        .def( "drawR", &FirstPassageGreensFunction::drawR )
        .def( "p_survival", &FirstPassageGreensFunction::p_survival )
        .def( "p_int_r", &FirstPassageGreensFunction::p_int_r )
        .def( "p_int_r_free", &FirstPassageGreensFunction::p_int_r_free )
        //.def( "p_r_fourier", &FirstPassageGreensFunction::p_r_fourier )
        .def( "dump", &FirstPassageGreensFunction::dump )
        ;



    class_<BasicPairGreensFunction>( "BasicPairGreensFunction",
                                     init<const Real, 
                                     const Real, 
                                     const Real>() )
        .def( "getD", &BasicPairGreensFunction::getD )
        .def( "getkf", &BasicPairGreensFunction::getkf )
        .def( "getSigma", &BasicPairGreensFunction::getSigma )
        .def( "drawTime", &BasicPairGreensFunction::drawTime )
        .def( "drawR", &BasicPairGreensFunction::drawR )
        .def( "drawTheta", &BasicPairGreensFunction::drawTheta )

//        .def( "p_tot", &BasicPairGreensFunction::p_tot )
        .def( "p_free", &BasicPairGreensFunction::p_free )
        .def( "ip_free", &BasicPairGreensFunction::ip_free )
        .def( "p_corr", &BasicPairGreensFunction::p_corr )
        .def( "ip_corr", &BasicPairGreensFunction::ip_corr )
        .def( "p_survival", &BasicPairGreensFunction::p_survival )
        .def( "p_int_r", &BasicPairGreensFunction::p_int_r )
        .def( "p_theta", &BasicPairGreensFunction::p_theta )
        .def( "ip_theta", &BasicPairGreensFunction::ip_theta )

        .def( "dump", &BasicPairGreensFunction::dump )
        ;

    class_<FreePairGreensFunction>( "FreePairGreensFunction",
                                    init<const Real>() )
        .def( "getD", &FreePairGreensFunction::getD )
        .def( "getkf", &FreePairGreensFunction::getkf )
        .def( "getSigma", &FreePairGreensFunction::getSigma )
        .def( "drawTime", &FreePairGreensFunction::drawTime )
        .def( "drawR", &FreePairGreensFunction::drawR )
        .def( "drawTheta", &FreePairGreensFunction::drawTheta )

        .def( "p_r", &FreePairGreensFunction::p_r )
        .def( "ip_r", &FreePairGreensFunction::ip_r )
        .def( "p_theta", &FreePairGreensFunction::p_theta )
        .def( "ip_theta", &FreePairGreensFunction::ip_theta )

        .def( "dump", &FreePairGreensFunction::dump )
        ;

    enum_<EventType>( "EventType" )
        .value( "REACTION", REACTION )
        .value( "ESCAPE", ESCAPE )
        ;

    class_<FirstPassagePairGreensFunction>( "FirstPassagePairGreensFunction",
                                            init<const Real, 
                                            const Real,
                                            const Real>() )
        .def( "seta", &FirstPassagePairGreensFunction::seta )
        .def( "geta", &FirstPassagePairGreensFunction::geta )
        .def( "getD", &FirstPassagePairGreensFunction::getD )
        .def( "getkf", &BasicPairGreensFunction::getkf )
        .def( "getSigma", &BasicPairGreensFunction::getSigma )
        .def( "drawTime", &FirstPassagePairGreensFunction::drawTime )
        //.def( "drawTime2", &FirstPassagePairGreensFunction::drawTime2 )
        .def( "drawEventType", &FirstPassagePairGreensFunction::drawEventType )
        .def( "drawR", &FirstPassagePairGreensFunction::drawR )
        .def( "drawTheta", &FirstPassagePairGreensFunction::drawTheta )

        .def( "p_survival", &FirstPassagePairGreensFunction::p_survival )
        .def( "dp_survival", &FirstPassagePairGreensFunction::dp_survival )
        .def( "p_leaves", &FirstPassagePairGreensFunction::p_leaves )
        .def( "p_leavea", &FirstPassagePairGreensFunction::p_leavea )
        .def( "leaves", &FirstPassagePairGreensFunction::leaves )
        .def( "leavea", &FirstPassagePairGreensFunction::leavea )

        .def( "p_0", &FirstPassagePairGreensFunction::p_0 )
        .def( "p_int_r", &FirstPassagePairGreensFunction::p_int_r )
        .def( "p_int_r", &FirstPassagePairGreensFunction::p_int_r )
        .def( "p_theta", &FirstPassagePairGreensFunction::p_theta )
        .def( "ip_theta", &FirstPassagePairGreensFunction::ip_theta )
        .def( "idp_theta", &FirstPassagePairGreensFunction::idp_theta )

        .def( "f_alpha0", &FirstPassagePairGreensFunction::f_alpha0 )
        .def( "alpha0_i", &FirstPassagePairGreensFunction::alpha0_i )
        .def( "f_alpha", &FirstPassagePairGreensFunction::f_alpha )
        .def( "f_alpha_aux", &FirstPassagePairGreensFunction::f_alpha_aux )

        .def( "p_survival_i_exp", &FirstPassagePairGreensFunction::p_survival_i_exp )
        .def( "p_survival_i_alpha", &FirstPassagePairGreensFunction::p_survival_i_alpha )

        //.def( "guess_maxi", &FirstPassagePairGreensFunction::guess_maxi )

        .def( "dump", &FirstPassagePairGreensFunction::dump )

//        .def( "alpha_i", &FirstPassagePairGreensFunction::alpha_i )
        ;


    class_<FirstPassageNoCollisionPairGreensFunction>
        ( "FirstPassageNoCollisionPairGreensFunction", init<const Real>() ) 
        .def( "seta", &FirstPassageNoCollisionPairGreensFunction::seta )
        .def( "geta", &FirstPassageNoCollisionPairGreensFunction::geta )
        .def( "getD", &FirstPassageNoCollisionPairGreensFunction::getD )
        .def( "drawTime", &FirstPassageNoCollisionPairGreensFunction::drawTime )
        .def( "drawR", &FirstPassageNoCollisionPairGreensFunction::drawR )
        .def( "drawTheta", 
              &FirstPassageNoCollisionPairGreensFunction::drawTheta )

        .def( "p_survival", 
              &FirstPassageNoCollisionPairGreensFunction::p_survival )
        .def( "dp_survival", 
              &FirstPassageNoCollisionPairGreensFunction::dp_survival )
        .def( "p_int_r", &FirstPassageNoCollisionPairGreensFunction::p_int_r )
        .def( "p_theta", &FirstPassageNoCollisionPairGreensFunction::p_theta )
        .def( "ip_theta", &FirstPassageNoCollisionPairGreensFunction::ip_theta )
        .def( "idp_theta", 
              &FirstPassageNoCollisionPairGreensFunction::idp_theta )

        .def( "dump", &FirstPassageNoCollisionPairGreensFunction::dump )
        ;

    def( "lengthSq", &length_sq< vector_type > );
    def( "length", &length< vector_type > );
    def( "distanceSq", &distance_sq< vector_type, vector_type > );
    def( "distance", &distance< vector_type, vector_type > );
    def( "distanceSq_Cyclic", &distance_sq_cyclic< vector_type, vector_type > );
    def( "distance_Cyclic", &distance_cyclic< vector_type, vector_type > );

    to_python_converter<position_type,
        position_to_ndarray_converter>();
    to_python_converter<Sphere<length_type>,
        sphere_to_python_converter>();
    peer::util::to_native_converter<Sphere<length_type>,
        python_to_sphere_converter>();
    peer::util::to_native_converter<position_type,
        ndarray_to_position_converter>();
    peer::util::to_native_converter<position_type,
        seq_to_position_converter>();

#if OBJECTMATRIX_USE_ITERATOR
    peer::util::register_stop_iteration_exc_translator();
#endif
    peer::MatrixSpace::__register_class();
    peer::SpeciesType::__register_class();

    class_<Model, boost::noncopyable>("Model")
        .add_property("network_rules",
            make_function(&Model::network_rules, return_internal_reference<>()))
        .def("new_species_type", &Model::new_species_type,
                return_internal_reference<>())
        ;

    peer::ReactionRule::__register_class();

    peer::IdentifierWrapper<SpeciesTypeID>::__register_class("SpeciesTypeID");
    peer::util::GeneratorIteratorWrapper<ptr_generator<NetworkRules::reaction_rule_generator> >::__register_class("reaction_rule_geneator");

    class_<NetworkRules, boost::noncopyable>("NetworkRules", no_init)
        .def("add_reaction_rule", &NetworkRules::add_reaction_rule)
        .def("query_reaction_rule", static_cast<NetworkRules::reaction_rule_generator*(NetworkRules::*)(SpeciesType const*)>(&NetworkRules::query_reaction_rule), return_value_policy<return_by_value>())
        .def("query_reaction_rule", static_cast<NetworkRules::reaction_rule_generator*(NetworkRules::*)(SpeciesType const*, SpeciesType const*)>(&NetworkRules::query_reaction_rule), return_value_policy<return_by_value>())
        ;
}
