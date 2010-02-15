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
#include <boost/python/reference_existing_object.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <numpy/arrayobject.h>

#include "PyEventScheduler.hpp"
#include "freeFunctions.hpp"
#include "FreeGreensFunction.hpp"
#include "FirstPassageGreensFunction.hpp"
#include "BasicPairGreensFunction.hpp"
#include "FreePairGreensFunction.hpp"
#include "FirstPassagePairGreensFunction.hpp"
#include "FirstPassageNoCollisionPairGreensFunction.hpp"

#include "utils/array_traits.hpp"
#include "MatrixSpace.hpp"
#include "Vector3.hpp"
#include "Sphere.hpp"
#include "utils.hpp"
#include "Model.hpp"
#include "World.hpp"
#include "EGFRDSimulator.hpp"
#include "NetworkRulesWrapper.hpp"
#include "ReactionRuleInfo.hpp"

#include "peer/utils.hpp"
#include "peer/py_hash_support.hpp"
#include "peer/tuple_converters.hpp"
#include "peer/set_indexing_suite.hpp"
#include "peer/numpy/wrapped_multi_array.hpp"
#include "peer/numpy/scalar_converters.hpp"
#include "peer/Particle.hpp"
#include "peer/SphericalShellWrapper.hpp"
#include "peer/CylindricalShellWrapper.hpp"
#include "peer/MatrixSpace.hpp"
#include "peer/SpeciesType.hpp"
#include "peer/Identifier.hpp"
#include "peer/ReactionRule.hpp"
#include "peer/GeneratorIteratorWrapper.hpp"
#include "peer/Exception.hpp"
#include "peer/GSLRandomNumberGenerator.hpp"

typedef CyclicWorldTraits<Real, Real> world_traits_type;
typedef World<world_traits_type> CyclicWorld;
typedef CyclicWorld::transaction_type transaction_type;
typedef EGFRDSimulatorTraitsBase<CyclicWorld> egfrd_simulator_traits_type;

static boost::python::object species_type_class;

struct position_to_ndarray_converter
{
    typedef world_traits_type::position_type native_type;
    
    static PyObject* convert( const native_type& p )
    {
        static const npy_intp dims[1] = { native_type::size() };
        void* data( PyDataMem_NEW( native_type::size() * sizeof( native_type::value_type ) ) );
        memcpy( data, static_cast<const void*>( p.data() ),
                native_type::size() * sizeof( native_type::value_type ) );
        PyObject* array( PyArray_New( &PyArray_Type, 1, 
                                      const_cast<npy_intp*>( dims ),
                                      peer::util::get_numpy_typecode<
                                          native_type::value_type >::value,
                                      NULL, data, 0, NPY_CARRAY, NULL ) );
        reinterpret_cast<PyArrayObject*>( array )->flags |= NPY_OWNDATA;
        return array;
    }
};

struct ndarray_to_position_converter
{
    typedef world_traits_type::position_type native_type;
    
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
                                     native_type::value_type >::value ), 0) );
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
            reinterpret_cast<native_type::value_type*>(PyArray_DATA(array_obj)));
        boost::python::decref(reinterpret_cast<PyObject*>(array_obj));
    }
};

struct seq_to_position_converter
{
    typedef world_traits_type::position_type native_type;
    
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
            PyFloat_AsDouble(boost::python::handle<>(PySequence_GetItem(ptr, 0)).get()),
            PyFloat_AsDouble(boost::python::handle<>(PySequence_GetItem(ptr, 1)).get()),
            PyFloat_AsDouble(boost::python::handle<>(PySequence_GetItem(ptr, 2)).get()));
    }
};

struct sphere_to_python_converter
{
    typedef Sphere<world_traits_type::length_type> native_type;

    static PyObject* convert(native_type const& v)
    {
        return boost::python::incref(
            boost::python::object(boost::make_tuple(
                v.position(), v.radius())).ptr());
    }
};

struct python_to_sphere_converter
{
    typedef Sphere<world_traits_type::length_type> native_type;

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

static Model::species_type_iterator Model_get_species_types_begin(Model& self)
{
    return self.get_species_types().begin();
}

static Model::species_type_iterator Model_get_species_types_end(Model& self)
{
    return self.get_species_types().end();
}

struct species_type_to_species_type_id_converter
{
    typedef world_traits_type::species_type_id_type native_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PyObject_TypeCheck(pyo, reinterpret_cast<PyTypeObject*>(
                species_type_class.ptr())))
        {
            return 0;
        }
        return pyo;
    }

    static void construct(PyObject* pyo, 
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        using namespace boost::python;
        void* storage(reinterpret_cast<
            converter::rvalue_from_python_storage<native_type>* >(
                data)->storage.bytes);
        new (storage) native_type(static_cast<SpeciesType*>(extract<SpeciesType*>(object(borrowed(pyo))))->id());
        data->convertible = storage;
    }
};

struct species_info_to_species_id_converter
{
    typedef world_traits_type::species_type::identifier_type native_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PyObject_TypeCheck(pyo, reinterpret_cast<PyTypeObject*>(
                species_type_class.ptr())))
        {
            return 0;
        }
        return pyo;
    }

    static void construct(PyObject* pyo, 
                          boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        using namespace boost::python;
        void* storage(reinterpret_cast<
            converter::rvalue_from_python_storage<native_type>* >(
                data)->storage.bytes);
        new (storage) native_type(static_cast<world_traits_type::species_type*>(extract<world_traits_type::species_type*>(object(borrowed(pyo))))->id());
        data->convertible = storage;
    }
};

static void Model___setitem__(Model* model, std::string const& key, std::string const& value)
{
    (*model)[key] = value;
}

template<typename T_>
struct world_get_species
{
    typedef typename T_::species_iterator iterator;

    static iterator begin(T_& impl)
    {
        return impl.get_species().begin();
    }

    static iterator end(T_& impl)
    {
        return impl.get_species().end();
    }
};

template<typename Tworld_>
struct SimulationStateTraits
{
    typedef Tworld_ world_type;
    typedef NetworkRulesWrapper<
            NetworkRules,
            ReactionRuleInfo<
                typename ReactionRule::identifier_type,
                SpeciesID,
                Real> > network_rules_type;
};


template<typename Ttraits_>
struct SimulationState
{
    Real get_dt() const
    {
        return boost::python::extract<Real>(self_.attr("dt"));
    }

    typename boost::shared_ptr<gsl_rng> get_rng()
    {
        return rng_;
    }

    typename Ttraits_::network_rules_type const& get_network_rules() const
    {
        return boost::python::extract<typename Ttraits_::network_rules_type const&>(self_.attr("network_rules"));
    }

    SimulationState(PyObject* self)
        : self_(boost::python::borrowed(self)),
          rng_(boost::shared_ptr<gsl_rng>(gsl_rng_alloc(gsl_rng_mt19937),
                                          gsl_rng_free))
        {}

private:
    boost::python::object self_;
    boost::shared_ptr<gsl_rng> rng_;
};


struct gsl_rng_to_python_converter
{
    typedef boost::shared_ptr<gsl_rng> native_type;

    static PyObject* convert(native_type const& v)
    {
        return boost::python::incref(
            boost::python::object(peer::GSLRandomNumberGenerator(v)).ptr());
    }
};


namespace boost { namespace python {

template<typename Ttraits_>
struct has_back_reference<SimulationState<Ttraits_> >: boost::mpl::true_ {};

} } // namespace boost::python

template<typename Twrapper_>
struct MatrixSpace_to_select_second_range_converter
{
    typedef Twrapper_ wrapper_type;
    typedef typename Twrapper_::impl_type impl_type;
    typedef select_first_range<const impl_type> native_type;

    static void* convertible(PyObject* ptr)
    {
        if (!PyObject_TypeCheck(ptr, wrapper_type::__class__))
        {
            return NULL;
        }

        return ptr;
    }
    
    static void construct(PyObject* ptr,
                          boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        Twrapper_& wrapper(boost::python::extract<Twrapper_&>(
            static_cast<PyObject*>(data->stage1.convertible)));
        data->stage1.convertible = new(data->storage.bytes) native_type(
            static_cast<impl_type&>(wrapper));
    }

    static void __register()
    {
        peer::util::to_native_converter<native_type,
            MatrixSpace_to_select_second_range_converter>();
    }
};


template<typename T_>
void register_id_generator(char const* class_name)
{
    using namespace boost::python;

    class_<SerialIDGenerator<T_> >(class_name, init<int>())
        .def("__call__", &SerialIDGenerator<T_>::operator())
        ;
}

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

    def( "normalize", &normalize<vector_type> );
    def( "cyclic_transpose", &cyclic_transpose<vector_type> );
    def( "calculate_pair_CoM", &calculate_pair_CoM<vector_type> );
    def( "apply_boundary", &apply_boundary<vector_type> );

    to_python_converter<world_traits_type::position_type,
        position_to_ndarray_converter>();
    to_python_converter<Sphere<world_traits_type::length_type>,
        sphere_to_python_converter>();
    peer::util::to_native_converter<Sphere<world_traits_type::length_type>,
        python_to_sphere_converter>();
    peer::util::to_native_converter<world_traits_type::position_type,
        ndarray_to_position_converter>();
    peer::util::to_native_converter<world_traits_type::position_type,
        seq_to_position_converter>();

    peer::MatrixSpace< MatrixSpace<egfrd_simulator_traits_type::spherical_shell_type, egfrd_simulator_traits_type::shell_id_type> >::__register_class("SphericalShellContainer");
    peer::MatrixSpace< MatrixSpace<egfrd_simulator_traits_type::cylindrical_shell_type, egfrd_simulator_traits_type::shell_id_type> >::__register_class("CylindricalShellContainer");
    peer::MatrixSpace<MatrixSpace<
        world_traits_type::particle_type,
        world_traits_type::particle_id_type, get_mapper_mf> >::__register_class("ParticleContainer");
    species_type_class = peer::SpeciesType::__register_class();

    class_<std::set<world_traits_type::particle_id_type> >("ParticleIDSet")
        .def(peer::util::set_indexing_suite<std::set<world_traits_type::particle_id_type> >())
        ;

    class_<Model, boost::noncopyable>("Model")
        .add_property("network_rules",
            make_function(&Model::network_rules,
                return_value_policy<reference_existing_object>()))
        .def("new_species_type", &Model::new_species_type,
                return_value_policy<reference_existing_object>())
        .def("get_species_type_by_id", &Model::get_species_type_by_id,
                return_value_policy<reference_existing_object>())
        .def("__getitem__", (std::string const&(Model::*)(std::string const&) const)
                &Model::operator[], return_value_policy<copy_const_reference>())
        .def("__setitem__", &Model___setitem__)
        .add_property("attributes",
                peer::util::range_from_range<
                    Model::attributes_range, Model, &Model::attributes>())
        .add_property("species_types",
                peer::util::range_from_range<
                    Model::species_type_range,
                    Model, &Model::get_species_types,
                    return_value_policy<reference_existing_object> >());
        ;

    peer::ReactionRule::__register_class();

    peer::IdentifierWrapper<world_traits_type::species_type_id_type>::__register_class("SpeciesTypeID");
    peer::util::to_native_converter<world_traits_type::species_type_id_type, species_type_to_species_type_id_converter>();

    peer::IdentifierWrapper<world_traits_type::species_id_type>::__register_class("SpeciesID");

    peer::util::GeneratorIteratorWrapper<ptr_generator<NetworkRules::reaction_rule_generator> >::__register_class("ReactionRuleGenerator");

    peer::util::ExceptionWrapper<not_found, peer::util::PyExcTraits<&PyExc_LookupError> >::__register_class("NotFound");
    peer::util::ExceptionWrapper<already_exists, peer::util::PyExcTraits<&PyExc_StandardError> >::__register_class("AlreadyExists");

    peer::IdentifierWrapper<world_traits_type::particle_id_type>::__register_class("ParticleID");
    register_id_generator<world_traits_type::particle_id_type>("ParticleIDGenerator");
    peer::ParticleWrapper<world_traits_type::particle_type>::__register_class("Particle");

    peer::IdentifierWrapper<egfrd_simulator_traits_type::shell_id_type>::__register_class("ShellID");
    register_id_generator<egfrd_simulator_traits_type::shell_id_type>("ShellIDGenerator");
    peer::SphericalShellWrapper<egfrd_simulator_traits_type::spherical_shell_type>::__register_class("SphericalShell");
    peer::CylindricalShellWrapper<egfrd_simulator_traits_type::cylindrical_shell_type>::__register_class("CylindricalShell");

    peer::IdentifierWrapper<egfrd_simulator_traits_type::domain_id_type>::__register_class("DomainID");
    register_id_generator<egfrd_simulator_traits_type::domain_id_type>("DomainIDGenerator");

    class_<NetworkRules, boost::noncopyable>("NetworkRules", no_init)
        .def("add_reaction_rule", &NetworkRules::add_reaction_rule)
        .def("query_reaction_rule", static_cast<NetworkRules::reaction_rule_generator*(NetworkRules::*)(SpeciesTypeID const&) const>(&NetworkRules::query_reaction_rule), return_value_policy<return_by_value>())
        .def("query_reaction_rule", static_cast<NetworkRules::reaction_rule_generator*(NetworkRules::*)(SpeciesTypeID const&, SpeciesTypeID const&) const>(&NetworkRules::query_reaction_rule), return_value_policy<return_by_value>())
        ;

    peer::util::register_tuple_converter<CyclicWorld::particle_id_pair>();

    peer::util::GeneratorIteratorWrapper<ptr_generator<CyclicWorld::particle_id_pair_generator> >::__register_class("ParticleIDPairGenerator");

    class_<transaction_type, boost::noncopyable>("Transaction", no_init)
        .add_property("num_particles", &transaction_type::num_particles)
        .add_property("added_particles",
            make_function(&transaction_type::get_added_particles,
                return_value_policy<return_by_value>()))
        .add_property("removed_particles",
            make_function(&transaction_type::get_removed_particles,
                return_value_policy<return_by_value>()))
        .add_property("modified_particles",
            make_function(&transaction_type::get_modified_particles,
                return_value_policy<return_by_value>()))
        .def("new_particle", &transaction_type::new_particle)
        .def("update_particle", &transaction_type::update_particle)
        .def("remove_particle", &transaction_type::remove_particle)
        .def("get_particle", &transaction_type::get_particle)
        .def("check_overlap", (transaction_type::particle_id_pair_list*(transaction_type::*)(transaction_type::particle_id_pair const&) const)&transaction_type::check_overlap, return_value_policy<manage_new_object>())
        .def("check_overlap", (transaction_type::particle_id_pair_list*(transaction_type::*)(transaction_type::sphere_type const&, transaction_type::particle_id_type const&) const)&transaction_type::check_overlap, return_value_policy<manage_new_object>())
        .def("check_overlap", (transaction_type::particle_id_pair_list*(transaction_type::*)(transaction_type::sphere_type const&) const)&transaction_type::check_overlap, return_value_policy<manage_new_object>())
        .def("create_transaction", &transaction_type::create_transaction,
                return_value_policy<manage_new_object>())
        .def("rollback", &transaction_type::rollback)
        .def("__iter__", &transaction_type::get_particles,
                return_value_policy<return_by_value>())
        ;

    class_<CyclicWorld>("World", init<CyclicWorld::length_type,
                                  CyclicWorld::size_type>())
        .add_property("num_particles", &CyclicWorld::num_particles)
        .add_property("world_size", &CyclicWorld::world_size)
        .add_property("cell_size", &CyclicWorld::cell_size)
        .add_property("matrix_size", &CyclicWorld::matrix_size)
        .add_property("species",
            range<return_value_policy<return_by_value>, CyclicWorld const&>(
                &world_get_species<CyclicWorld const>::begin,
                &world_get_species<CyclicWorld const>::end))
        .def("add_species", &CyclicWorld::add_species)
        .def("get_species",
            (CyclicWorld::species_type const&(CyclicWorld::*)(CyclicWorld::species_id_type const&) const)&CyclicWorld::get_species,
            return_internal_reference<>())
        .def("distance", &CyclicWorld::distance)
        .def("distance_sq", &CyclicWorld::distance_sq)
        .def("apply_boundary", (CyclicWorld::position_type(CyclicWorld::*)(CyclicWorld::position_type const&) const)&CyclicWorld::apply_boundary)
        .def("apply_boundary", (CyclicWorld::length_type(CyclicWorld::*)(CyclicWorld::length_type const&) const)&CyclicWorld::apply_boundary)
        .def("cyclic_transpose", (CyclicWorld::position_type(CyclicWorld::*)(CyclicWorld::position_type const&, CyclicWorld::position_type const&) const)&CyclicWorld::cyclic_transpose)
        .def("cyclic_transpose", (CyclicWorld::length_type(CyclicWorld::*)(CyclicWorld::length_type const&, CyclicWorld::length_type const&) const)&CyclicWorld::cyclic_transpose)
        .def("new_particle", &CyclicWorld::new_particle)
        .def("check_overlap", (CyclicWorld::particle_id_pair_list*(CyclicWorld::*)(CyclicWorld::particle_id_pair const&) const)&CyclicWorld::check_overlap, return_value_policy<manage_new_object>())
        .def("check_overlap", (CyclicWorld::particle_id_pair_list*(CyclicWorld::*)(CyclicWorld::sphere_type const&, CyclicWorld::particle_id_type const&) const)&CyclicWorld::check_overlap, return_value_policy<manage_new_object>())
        .def("check_overlap", (CyclicWorld::particle_id_pair_list*(CyclicWorld::*)(CyclicWorld::sphere_type const&) const)&CyclicWorld::check_overlap, return_value_policy<manage_new_object>())
        .def("update_particle", &CyclicWorld::update_particle)
        .def("remove_particle", &CyclicWorld::remove_particle)
        .def("get_particle", &CyclicWorld::get_particle)
        .def("create_transaction", &CyclicWorld::create_transaction,
                return_value_policy<manage_new_object>())
        .def("__iter__", &CyclicWorld::get_particles,
                return_value_policy<return_by_value>())
        ;

    typedef world_traits_type::species_type species_type;

    class_<species_type>("SpeciesInfo",
            init<species_type::identifier_type, species_type::species_type_id_type>())
        .def(init<species_type::identifier_type, species_type::species_type_id_type, species_type::length_type, species_type::D_type>())
        .add_property("id",
            make_function(&species_type::id,
                return_value_policy<return_by_value>()))
        .add_property("type_id",
            make_function(&species_type::type_id,
                return_value_policy<return_by_value>()))
        .add_property("radius",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    species_type, length_type,
                    &species_type::radius,
                    &species_type::radius>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                species_type, Real,
                &species_type::radius,
                &species_type::radius>::set)
        ;

    peer::util::to_native_converter<world_traits_type::species_id_type, species_info_to_species_id_converter>();

    typedef SimulationStateTraits<CyclicWorld> simulation_state_traits_type;
    typedef SimulationState<simulation_state_traits_type> simulation_state_type;
    class_<simulation_state_type, boost::noncopyable>("SimulationState")
        .add_property("rng", &simulation_state_type::get_rng)
        ;

    peer::GSLRandomNumberGenerator::__register_class<gsl_rng_mt19937>("RandomNumberGenerator");

    to_python_converter<boost::shared_ptr<gsl_rng>,
            gsl_rng_to_python_converter>();

    peer::util::register_scalar_to_native_converters();
}
