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


#include <numpy/arrayobject.h>

#include "distance.hpp"

#include "wrapped_multi_array.hpp"

#include "PyEventScheduler.hpp"
#include "freeFunctions.hpp"
#include "FreeGreensFunction.hpp"
#include "FirstPassageGreensFunction.hpp"
#include "BasicPairGreensFunction.hpp"
#include "FreePairGreensFunction.hpp"
#include "FirstPassagePairGreensFunction.hpp"
#include "FirstPassageNoCollisionPairGreensFunction.hpp"



using namespace boost::python;

boost::python::tuple tuple_to_python( boost::tuples::null_type )
{
    return boost::python::tuple();
}

template <class H, class T>
boost::python::tuple tuple_to_python( const boost::tuples::cons<H,T>& x )
{
    return boost::python::tuple( boost::python::make_tuple( x.get_head() ) +
                                 tuple_to_python( x.get_tail() ) );
}

template <class T>
struct tupleconverter
{
    static PyObject* convert( const T& x )
    {
        return incref( tuple_to_python(x).ptr() );
    }
};

void* extract_pyarray(PyObject* x)
{
    return PyObject_TypeCheck(x, &PyArray_Type) ? x : 0;
}






template<typename T_, std::size_t Ndims_>
class ndarray_wrapped_multi_array_converter
{
public:
    typedef python_array_lifecycle_manager<T_> lcmgr_type;
    typedef wrapped_multi_array<T_, Ndims_> native_type;

public:
    static void* convertible(PyObject* ptr)
    {
        if (!PyArray_Check(ptr))
        {
            return NULL;
        }

        PyObject* retval(
            PyArray_CastToType(
                reinterpret_cast<PyArrayObject*>(ptr),
                PyArray_DescrFromType(
                    get_numpy_typecode<
                        typename native_type::element>::value), 0));
        if (!retval)
        {
            return NULL;
        }

        if (PyArray_NDIM(reinterpret_cast<PyArrayObject*>(retval)) != Ndims_)
        {
            boost::python::decref(retval);
            return NULL;
        }

#ifdef DEBUG
        if (retval != ptr)
        {
            std::cerr << "copy performed" << std::endl;
        }
#endif

        return retval;
    }

    static void construct(PyObject* ptr,
            boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        PyArrayObject* array_obj = static_cast<PyArrayObject*>(
                data->stage1.convertible);
        typename native_type::index_list ma_strides;

        for (std::size_t i = 0; i < Ndims_; ++i)
        {
            ma_strides[i] = array_obj->strides[i] / sizeof(T_);
        }

        data->stage1.convertible = new(data->storage.bytes) native_type(
                new lcmgr_type(array_obj),
                *reinterpret_cast<boost::array<npy_intp, Ndims_>*>(
                    array_obj->dimensions),
                ma_strides,
                PyArray_ISCONTIGUOUS(array_obj) ?
                    static_cast<boost::general_storage_order<Ndims_> >(
                        boost::c_storage_order()):
                    static_cast<boost::general_storage_order<Ndims_> >(
                        boost::fortran_storage_order())
                );
    }
};

template<typename T_, std::size_t Ndims_>
void register_ndarray_wrapped_multi_array_converter()
{
    typedef ndarray_wrapped_multi_array_converter<T_, Ndims_> Converter;
    boost::python::converter::registry::push_back(
        &Converter::convertible,
        reinterpret_cast<boost::python::converter::constructor_function>(
            &Converter::construct),
        boost::python::type_id<typename Converter::native_type>());
}




// Exception translators here.

void translateException( const std::exception& anException )
{
  PyErr_SetString( PyExc_RuntimeError, anException.what() );
}


// GSL error handler.
void __gsl_error_handler( const char* reason,
			     const char* file,
			     int line,
			     int gsl_errno )
{
    throw std::runtime_error( std::string( "GSL error: " ) +
			      std::string( reason ) +
			      std::string( " at " ) +
			      std::string( file ) + std::string( ":" ) +
			      boost::lexical_cast<std::string>( line ) );
}


BOOST_PYTHON_MODULE( _gfrd )
{
    import_array();

    // GSL error handler: is this the best place for this?
    gsl_set_error_handler( &__gsl_error_handler );

  
    register_exception_translator<std::exception>( &translateException );


    register_ndarray_wrapped_multi_array_converter<npy_double, 1>();
    register_ndarray_wrapped_multi_array_converter<npy_double, 2>();
    register_ndarray_wrapped_multi_array_converter<npy_double, 3>();


//    to_python_converter<PyEvent, PyEvent_to_python>();

    // boost::tuple
    to_python_converter<boost::tuple<Real,EventType>, 
        tupleconverter<boost::tuple<Real,EventType> > >();


    // free functions
    //def( "p_free", p_free );
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
//	.def( "fire", &PyEvent::fire )
//	.def( "update", &PyEvent::update )
//	.def( "isDependentOn", &PyEvent::isDependentOn )
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
//	.def( "updateAllEventDependency", 
//	      &PyEventScheduler::updateAllEventDependency )
	.def( "check", &PyEventScheduler::check )
	;


    //    converter::registry::insert( &extract_pyarray, type_id<PyArrayObject>());


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

//	.def( "alpha_i", &FirstPassagePairGreensFunction::alpha_i )
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

//	.def( "alpha_i", &FirstPassageNoCollisionPairGreensFunction::alpha_i )
	;

    def( "distanceSq", &distanceSq );
    def( "distance", &distance );
    def( "distanceSq_Cyclic", &distanceSq_Cyclic );
    def( "distance_Cyclic", &distance_Cyclic );

}
