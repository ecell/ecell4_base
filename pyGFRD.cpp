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
#include "utils.hpp"
#include "geometry.hpp"
#include "MatrixSpace.hpp"
#include "Vector3.hpp"
#include "Sphere.hpp"
#include "Cylinder.hpp"
#include "Point.hpp"
#include "Model.hpp"
#include "World.hpp"
#include "GSLRandomNumberGenerator.hpp"
#include "EGFRDSimulator.hpp"

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
#include "peer/RandomNumberGenerator.hpp"
#include "peer/STLIteratorWrapper.hpp"

typedef CyclicWorldTraits<Real, Real> world_traits_type;
typedef World<world_traits_type> CyclicWorld;
typedef Model model_type;
typedef CyclicWorld::transaction_type transaction_type;
typedef EGFRDSimulatorTraitsBase<CyclicWorld, model_type> egfrd_simulator_traits_type;

static boost::python::object species_type_class;
static boost::python::object species_info_class;

template<typename T_>
struct array_to_ndarray_converter
{
    typedef T_ native_type;
    
    static PyObject* convert( const native_type& p )
    {
        typedef typename boost::range_value<native_type>::type value_type;
        static const npy_intp dims[1] = { native_type::size() };
        void* data( PyDataMem_NEW( boost::size(p) * sizeof( value_type ) ) );
        memcpy( data, static_cast<const void*>( &p[0] ),
                native_type::size() * sizeof( value_type ) );
        PyObject* array( PyArray_New( &PyArray_Type, 1, 
                                      const_cast<npy_intp*>( dims ),
                                      peer::util::get_numpy_typecode<
                                          value_type >::value,
                                      NULL, data, 0, NPY_CARRAY, NULL ) );
        reinterpret_cast<PyArrayObject*>( array )->flags |= NPY_OWNDATA;
        return array;
    }
};

typedef array_to_ndarray_converter<world_traits_type::position_type> position_to_ndarray_converter;

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

struct species_type_to_species_id_converter
{
    typedef world_traits_type::species_id_type native_type;

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
                species_info_class.ptr())))
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

template<typename Ttraits_>
struct SimulationState
{
    typedef typename Ttraits_::time_type time_type;
    typedef typename Ttraits_::rng_type rng_type;
    typedef typename Ttraits_::network_rules_type network_rules_type;

    time_type get_dt() const
    {
        return boost::python::extract<Real>(self_.attr("dt"));
    }

    rng_type get_rng()
    {
        return rng_;
    }

    network_rules_type const& get_network_rules() const
    {
        return boost::python::extract<typename Ttraits_::network_rules_type const&>(self_.attr("network_rules"));
    }

    SimulationState(PyObject* self)
        : self_(boost::python::borrowed(self)) {}

private:
    boost::python::object self_;
    rng_type rng_;
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
static void register_id_generator(char const* class_name)
{
    using namespace boost::python;

    class_<SerialIDGenerator<T_> >(class_name, init<int>())
        .def("__call__", &SerialIDGenerator<T_>::operator())
        ;
}

template<typename T_>
static T_ 
calculate_pair_CoM(T_ const& p1, 
                   T_ const& p2, 
                   typename element_type_of< T_ >::type const& D1,
                   typename element_type_of< T_ >::type const& D2,
                   typename element_type_of< T_ >::type const& world_size)
{
    typedef typename element_type_of<T_>::type element_type;   

    T_ retval;

    const T_ p2t(cyclic_transpose<T_>(p2, p1, world_size));

    return modulo(
        divide(
            add(multiply(p1, D2), multiply(p2t, D1)),
            add(D1, D2)),
        world_size);
}

template<typename T_>
static boost::python::object Sphere___getitem__(T_ const& obj, int index)
{
    switch (index)
    {
    default:
        PyErr_SetString(PyExc_IndexError, "index out of range");
        boost::python::throw_error_already_set();
        break;
    case 0:
        return boost::python::object(obj.position());
    case 1:
        return boost::python::object(obj.radius());
    }

    return boost::python::object();
}

template<gsl_rng_type const*& Prng_>
static GSLRandomNumberGenerator create_gsl_rng()
{
    return GSLRandomNumberGenerator(gsl_rng_alloc(Prng_));
}

BOOST_PYTHON_MODULE( _gfrd )
{
    using namespace boost::python;

    import_array();

    // GSL error handler: is this the best place for this?
    gsl_set_error_handler( &gsl_error_handler );

  
    peer::util::register_std_exception_translator();

    peer::util::register_seq_wrapped_multi_array_converter<world_traits_type::length_type>();
    // peer::util::register_ndarray_wrapped_multi_array_converter<length_type, 1>();
    peer::util::register_ndarray_wrapped_multi_array_converter<world_traits_type::length_type, 2>();
    peer::util::register_ndarray_wrapped_multi_array_converter<world_traits_type::length_type, 3>();

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
        .value( "SINGLE_REACTION", SINGLE_REACTION )
        .value( "NOT_A_SINGLE_REACTION", NOT_A_SINGLE_REACTION )
        .value( "SINGLE_ESCAPE", SINGLE_ESCAPE )
        .value( "COM_ESCAPE", COM_ESCAPE )
        .value( "IV_ESCAPE", IV_ESCAPE )
        .value( "PAIR_REACTION", PAIR_REACTION )
        .value( "INTERACTION", INTERACTION )
        .value( "BURST", BURST )
        .value( "MULTI_ESCAPE", MULTI_ESCAPE )
        .value( "MULTI_REACTION", MULTI_REACTION )
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

    def( "length_sq", &length_sq< world_traits_type::position_type > );
    def( "length", &length< world_traits_type::position_type > );
    def( "distance", (world_traits_type::length_type(*)(world_traits_type::position_type const&, world_traits_type::position_type const&))&distance<world_traits_type::position_type> );
    def( "distance_cyclic", &distance_cyclic<world_traits_type::position_type, world_traits_type::position_type> );
    def( "apply_boundary", &apply_boundary<world_traits_type::position_type, world_traits_type::length_type> );
    def( "calculate_pair_CoM", &calculate_pair_CoM<world_traits_type::position_type> );

    def( "normalize", &normalize<world_traits_type::position_type> );
    def( "cyclic_transpose", &cyclic_transpose<world_traits_type::position_type, element_type_of<world_traits_type::position_type>::type> );

    to_python_converter<world_traits_type::position_type,
        position_to_ndarray_converter>();
#if 0
    to_python_converter<Sphere<world_traits_type::length_type>,
        sphere_to_python_converter>();
#endif
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

    class_<model_type, boost::noncopyable>("Model")
        .add_property("network_rules",
            make_function(&model_type::network_rules,
                return_value_policy<reference_existing_object>()))
        .def("new_species_type", &model_type::new_species_type,
                return_value_policy<reference_existing_object>())
        .def("get_species_type_by_id", &model_type::get_species_type_by_id,
                return_value_policy<reference_existing_object>())
        .def("__getitem__", (std::string const&(model_type::*)(std::string const&) const)
                &model_type::operator[], return_value_policy<copy_const_reference>())
        .def("__setitem__", &Model___setitem__)
        .add_property("attributes",
                peer::util::range_from_range<
                    model_type::attributes_range,
                    model_type, &model_type::attributes>())
        .add_property("species_types",
                peer::util::range_from_range<
                    model_type::species_type_range,
                    model_type, &model_type::get_species_types,
                    return_value_policy<reference_existing_object> >());
        ;

    peer::ReactionRule::__register_class();

    peer::IdentifierWrapper<world_traits_type::species_id_type>::__register_class("SpeciesTypeID");
    peer::util::to_native_converter<world_traits_type::species_id_type, species_type_to_species_id_converter>();

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
        .def("query_reaction_rule", static_cast<NetworkRules::reaction_rule_generator*(NetworkRules::*)(world_traits_type::species_id_type const&) const>(&NetworkRules::query_reaction_rule), return_value_policy<return_by_value>())
        .def("query_reaction_rule", static_cast<NetworkRules::reaction_rule_generator*(NetworkRules::*)(world_traits_type::species_id_type const&, world_traits_type::species_id_type const&) const>(&NetworkRules::query_reaction_rule), return_value_policy<return_by_value>())
        ;


    typedef egfrd_simulator_traits_type::reaction_rule_type reaction_rule_info_type;
    class_<reaction_rule_info_type>("ReactionRuleInfo", no_init)
        .add_property("id", 
            make_function(&reaction_rule_info_type::id,
                          return_value_policy<return_by_value>()))
        .add_property("k", make_function(&reaction_rule_info_type::k))
        .add_property("products",
            make_function(&reaction_rule_info_type::get_products,
                          return_value_policy<return_by_value>()))
        .add_property("reactants",
            make_function(&reaction_rule_info_type::get_reactants,
                          return_value_policy<return_by_value>()));

    peer::util::register_range_to_tuple_converter<reaction_rule_info_type::species_id_range>();

    typedef egfrd_simulator_traits_type::network_rules_type network_rules_wrapper_type;
    class_<network_rules_wrapper_type, boost::noncopyable>("NetworkRulesWrapper", init<NetworkRules const&>())
        .def("query_reaction_rule", (network_rules_wrapper_type::reaction_rule_vector const&(network_rules_wrapper_type::*)(network_rules_wrapper_type::species_id_type const&) const)&network_rules_wrapper_type::query_reaction_rule,
            return_value_policy<return_by_value>())
        .def("query_reaction_rule", (network_rules_wrapper_type::reaction_rule_vector const&(network_rules_wrapper_type::*)(network_rules_wrapper_type::species_id_type const&, network_rules_wrapper_type::species_id_type const&) const)&network_rules_wrapper_type::query_reaction_rule,
            return_value_policy<return_by_value>())
        ;
    peer::util::register_stl_iterator_range_converter<network_rules_wrapper_type::reaction_rule_vector>();

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
        .def("check_overlap", (transaction_type::particle_id_pair_list*(transaction_type::*)(transaction_type::particle_type::shape_type const&, transaction_type::particle_id_type const&) const)&transaction_type::check_overlap, return_value_policy<manage_new_object>())
        .def("check_overlap", (transaction_type::particle_id_pair_list*(transaction_type::*)(transaction_type::particle_type::shape_type const&) const)&transaction_type::check_overlap, return_value_policy<manage_new_object>())
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
            range<return_value_policy<return_by_value>, CyclicWorld const>(
                &world_get_species<CyclicWorld const>::begin,
                &world_get_species<CyclicWorld const>::end))
        .def("add_species", &CyclicWorld::add_species)
        .def("get_species",
            (CyclicWorld::species_type const&(CyclicWorld::*)(CyclicWorld::species_id_type const&) const)&CyclicWorld::get_species,
            return_internal_reference<>())
        .def("distance", &CyclicWorld::distance<egfrd_simulator_traits_type::sphere_type>)
        .def("distance", &CyclicWorld::distance<egfrd_simulator_traits_type::cylinder_type>)
        .def("distance", &CyclicWorld::distance<egfrd_simulator_traits_type::box_type>)
        .def("distance", &CyclicWorld::distance<CyclicWorld::position_type>)
        .def("apply_boundary", (CyclicWorld::position_type(CyclicWorld::*)(CyclicWorld::position_type const&) const)&CyclicWorld::apply_boundary)
        .def("apply_boundary", (CyclicWorld::length_type(CyclicWorld::*)(CyclicWorld::length_type const&) const)&CyclicWorld::apply_boundary)
        .def("cyclic_transpose", (CyclicWorld::position_type(CyclicWorld::*)(CyclicWorld::position_type const&, CyclicWorld::position_type const&) const)&CyclicWorld::cyclic_transpose)
        .def("cyclic_transpose", (CyclicWorld::length_type(CyclicWorld::*)(CyclicWorld::length_type const&, CyclicWorld::length_type const&) const)&CyclicWorld::cyclic_transpose)
        .def("calculate_pair_CoM", &CyclicWorld::calculate_pair_CoM<CyclicWorld::position_type>)
        .def("new_particle", &CyclicWorld::new_particle)
        .def("check_overlap", (CyclicWorld::particle_id_pair_list*(CyclicWorld::*)(CyclicWorld::particle_id_pair const&) const)&CyclicWorld::check_overlap, return_value_policy<manage_new_object>())
        .def("check_overlap", (CyclicWorld::particle_id_pair_list*(CyclicWorld::*)(CyclicWorld::particle_shape_type const&, CyclicWorld::particle_id_type const&) const)&CyclicWorld::check_overlap, return_value_policy<manage_new_object>())
        .def("check_overlap", (CyclicWorld::particle_id_pair_list*(CyclicWorld::*)(CyclicWorld::particle_shape_type const&) const)&CyclicWorld::check_overlap, return_value_policy<manage_new_object>())
        .def("update_particle", &CyclicWorld::update_particle)
        .def("remove_particle", &CyclicWorld::remove_particle)
        .def("get_particle", &CyclicWorld::get_particle)
        .def("create_transaction", &CyclicWorld::create_transaction,
                return_value_policy<manage_new_object>())
        .def("__iter__", &CyclicWorld::get_particles,
                return_value_policy<return_by_value>())
        ;

    typedef world_traits_type::species_type species_type;

    species_info_class = class_<species_type>("SpeciesInfo",
            init<species_type::identifier_type>())
        .def(init<species_type::identifier_type, species_type::length_type, species_type::D_type, species_type::surface_type>())
        .add_property("id",
            make_function(&species_type::id,
                return_value_policy<return_by_value>()))
        .add_property("radius",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    species_type, species_type::length_type,
                    &species_type::radius,
                    &species_type::radius>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                species_type, species_type::length_type,
                &species_type::radius,
                &species_type::radius>::set)
        .add_property("surface",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    species_type, species_type::surface_type,
                    &species_type::surface,
                    &species_type::surface>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                species_type, species_type::surface_type,
                &species_type::surface,
                &species_type::surface>::set)
        .add_property("D",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    species_type, species_type::D_type,
                    &species_type::D,
                    &species_type::D>::get,
                return_value_policy<return_by_value>()),
            &peer::util::reference_accessor_wrapper<
                species_type, species_type::D_type,
                &species_type::D,
                &species_type::D>::set)
        ;

    class_<egfrd_simulator_traits_type::sphere_type>("Sphere")
        .add_property("position",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::sphere_type,
                    egfrd_simulator_traits_type::sphere_type::position_type,
                    &egfrd_simulator_traits_type::sphere_type::position,
                    &egfrd_simulator_traits_type::sphere_type::position>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::sphere_type,
                    egfrd_simulator_traits_type::sphere_type::position_type,
                    &egfrd_simulator_traits_type::sphere_type::position,
                    &egfrd_simulator_traits_type::sphere_type::position>::set))
        .add_property("radius",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::sphere_type,
                    egfrd_simulator_traits_type::sphere_type::length_type,
                    &egfrd_simulator_traits_type::sphere_type::radius,
                    &egfrd_simulator_traits_type::sphere_type::radius>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::sphere_type,
                    egfrd_simulator_traits_type::sphere_type::length_type,
                    &egfrd_simulator_traits_type::sphere_type::radius,
                    &egfrd_simulator_traits_type::sphere_type::radius>::set))
        .def("__getitem__", &Sphere___getitem__<egfrd_simulator_traits_type::sphere_type>);

    class_<egfrd_simulator_traits_type::cylinder_type>("Cylinder")
        .add_property("position",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::position_type,
                    &egfrd_simulator_traits_type::cylinder_type::position,
                    &egfrd_simulator_traits_type::cylinder_type::position>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::position_type,
                    &egfrd_simulator_traits_type::cylinder_type::position,
                    &egfrd_simulator_traits_type::cylinder_type::position>::set))
        .add_property("radius",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::length_type,
                    &egfrd_simulator_traits_type::cylinder_type::radius,
                    &egfrd_simulator_traits_type::cylinder_type::radius>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::length_type,
                    &egfrd_simulator_traits_type::cylinder_type::radius,
                    &egfrd_simulator_traits_type::cylinder_type::radius>::set))
        .add_property("size",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::length_type,
                    &egfrd_simulator_traits_type::cylinder_type::size,
                    &egfrd_simulator_traits_type::cylinder_type::size>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::length_type,
                    &egfrd_simulator_traits_type::cylinder_type::size,
                    &egfrd_simulator_traits_type::cylinder_type::size>::set))
        .add_property("unit_z",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::position_type,
                    &egfrd_simulator_traits_type::cylinder_type::unit_z,
                    &egfrd_simulator_traits_type::cylinder_type::unit_z>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::cylinder_type,
                    egfrd_simulator_traits_type::cylinder_type::position_type,
                    &egfrd_simulator_traits_type::cylinder_type::unit_z,
                    &egfrd_simulator_traits_type::cylinder_type::unit_z>::set));

    to_python_converter<boost::array<egfrd_simulator_traits_type::box_type::length_type, 3>, array_to_ndarray_converter<boost::array<egfrd_simulator_traits_type::box_type::length_type, 3> > >();

    class_<egfrd_simulator_traits_type::box_type>("Box")
        .def(init<egfrd_simulator_traits_type::box_type::position_type, 
                  egfrd_simulator_traits_type::box_type::position_type, 
                  egfrd_simulator_traits_type::box_type::position_type, 
                  egfrd_simulator_traits_type::box_type::position_type, 
                  egfrd_simulator_traits_type::box_type::length_type,
                  egfrd_simulator_traits_type::box_type::length_type, 
                  egfrd_simulator_traits_type::box_type::length_type>())
        .add_property("position",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::position,
                    &egfrd_simulator_traits_type::box_type::position>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::position,
                    &egfrd_simulator_traits_type::box_type::position>::set))
        .add_property("extent",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    boost::array<egfrd_simulator_traits_type::box_type::length_type, 3>,
                    &egfrd_simulator_traits_type::box_type::extent,
                    &egfrd_simulator_traits_type::box_type::extent>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    boost::array<egfrd_simulator_traits_type::box_type::length_type, 3>,
                    &egfrd_simulator_traits_type::box_type::extent,
                    &egfrd_simulator_traits_type::box_type::extent>::set))
        .add_property("unit_x",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::unit_x,
                    &egfrd_simulator_traits_type::box_type::unit_x>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::unit_x,
                    &egfrd_simulator_traits_type::box_type::unit_x>::set))
        .add_property("unit_y",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::unit_y,
                    &egfrd_simulator_traits_type::box_type::unit_y>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::unit_y,
                    &egfrd_simulator_traits_type::box_type::unit_y>::set))
        .add_property("unit_z",
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::unit_z,
                    &egfrd_simulator_traits_type::box_type::unit_z>::get,
                return_value_policy<return_by_value>()),
            make_function(
                &peer::util::reference_accessor_wrapper<
                    egfrd_simulator_traits_type::box_type,
                    egfrd_simulator_traits_type::box_type::position_type,
                    &egfrd_simulator_traits_type::box_type::unit_z,
                    &egfrd_simulator_traits_type::box_type::unit_z>::set));

    peer::util::to_native_converter<world_traits_type::species_id_type, species_info_to_species_id_converter>();

    typedef SimulationState<egfrd_simulator_traits_type> simulation_state_type;
    class_<simulation_state_type, boost::noncopyable>("SimulationState")
        .add_property("rng", &simulation_state_type::get_rng)
        ;

    peer::RandomNumberGenerator<GSLRandomNumberGenerator>::__register_class("RandomNumberGenerator");
    def("create_gsl_rng", &create_gsl_rng<gsl_rng_mt19937>);

    peer::util::register_scalar_to_native_converters();
}
