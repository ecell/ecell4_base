#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/python.hpp>

#include "freeFunctions.hpp"
#include "FreeGreensFunction.hpp"
#include "FirstPassageGreensFunction.hpp"
#include "BasicPairGreensFunction.hpp"
#include "FreePairGreensFunction.hpp"
#include "FirstPassagePairGreensFunction.hpp"
#include "FirstPassageNoCollisionPairGreensFunction.hpp"

BOOST_PYTHON_MODULE( _greens_functions )
{
    using namespace boost::python;

    //import_array();

    enum_<EventType>( "EventType" )
        .value( "SINGLE_REACTION", SINGLE_REACTION )
        .value( "SINGLE_ESCAPE", SINGLE_ESCAPE )
        .value( "COM_ESCAPE", COM_ESCAPE )
        .value( "IV_EVENT", IV_EVENT )
        .value( "IV_ESCAPE", IV_ESCAPE )
        .value( "IV_REACTION", IV_REACTION )
        .value( "IV_INTERACTION", IV_INTERACTION )
        .value( "BURST", BURST )
        .value( "MULTI_ESCAPE", MULTI_ESCAPE )
        .value( "MULTI_REACTION", MULTI_REACTION )
        ;

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

    class_<FreeGreensFunction>("FreeGreensFunction", init<Real>())
        .def( "getD", &FreeGreensFunction::getD )
        .def( "drawTime", &FreeGreensFunction::drawTime )
        .def( "drawR", &FreeGreensFunction::drawR )
        .def( "p_r", &FreeGreensFunction::p_r )
        .def( "ip_r", &FreeGreensFunction::ip_r )
        .def( "dump", &FreeGreensFunction::dump )
        ;

    class_<FirstPassageGreensFunction>("FirstPassageGreensFunction",
                                       init<Real, Real>())
        .def( "getD", &FirstPassageGreensFunction::getD )
        .def( "geta", &FirstPassageGreensFunction::geta )
        .def( "drawTime", &FirstPassageGreensFunction::drawTime )
        .def( "drawR", &FirstPassageGreensFunction::drawR )
        .def( "p_survival", &FirstPassageGreensFunction::p_survival )
        .def( "p_int_r", &FirstPassageGreensFunction::p_int_r )
        .def( "p_int_r_free", &FirstPassageGreensFunction::p_int_r_free )
        //.def( "p_r_fourier", &FirstPassageGreensFunction::p_r_fourier )
        .def( "dump", &FirstPassageGreensFunction::dump )
        ;

    class_<BasicPairGreensFunction>("BasicPairGreensFunction",
                                    init<Real, Real, Real>())
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

    class_<FreePairGreensFunction>("FreePairGreensFunction", init<Real>())
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

    class_<FirstPassagePairGreensFunction>( "FirstPassagePairGreensFunction",
                                            init<Real, Real, Real, Real>() )
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
        ( "FirstPassageNoCollisionPairGreensFunction", init<Real, Real>() ) 
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
}
