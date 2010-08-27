#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/python.hpp>

#include "freeFunctions.hpp"
#include "GreensFunction1DAbsAbs.hpp"
#include "GreensFunction1DRadAbs.hpp"
#include "GreensFunction3DSym.hpp"
#include "GreensFunction3DAbsSym.hpp"
#include "GreensFunction3DRadInf.hpp"
#include "GreensFunction3D.hpp"
#include "GreensFunction3DRadAbs.hpp"
#include "GreensFunction3DAbs.hpp"

BOOST_PYTHON_MODULE( _greens_functions )
{
    using namespace boost::python;

    //import_array();
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

    
    class_<GreensFunction1DAbsAbs>( "GreensFunction1DAbsAbs",
					init<Real, Real, Real, Real>() )
	.def( init<Real, Real, Real, Real, Real>())
	.def( "getD", &FirstPassageGreensFunction1D::getD )
	.def( "getv", &FirstPassageGreensFunction1D::getv )
	.def( "getsigma", &FirstPassageGreensFunction1D::getsigma )
	.def( "seta", &FirstPassageGreensFunction1D::seta )
	.def( "geta", &FirstPassageGreensFunction1D::geta )
	.def( "setr0", &FirstPassageGreensFunction1D::setr0 )
	.def( "getr0", &FirstPassageGreensFunction1D::getr0 )
	.def( "drawTime", &FirstPassageGreensFunction1D::drawTime )
	.def( "drawR", &FirstPassageGreensFunction1D::drawR )
	.def( "drawEventType", &FirstPassageGreensFunction1D::drawEventType )
	.def( "leaves", &FirstPassageGreensFunction1D::leaves )
	.def( "leavea", &FirstPassageGreensFunction1D::leavea )
	.def( "p_survival", &FirstPassageGreensFunction1D::p_survival )
	.def( "calcpcum", &FirstPassageGreensFunction1D::calcpcum )
	;

    class_<GreensFunction1DRadAbs>( "GreensFunction1DRadAbs",
					init<Real, Real, Real, Real, Real>() )
	.def( init<Real, Real, Real, Real, Real, Real>())
	.def( "getk", &FirstPassageGreensFunction1DRad::getk )
	.def( "getD", &FirstPassageGreensFunction1DRad::getD )
	.def( "getv", &FirstPassageGreensFunction1DRad::getv )
	.def( "getsigma", &FirstPassageGreensFunction1DRad::getsigma )
	.def( "seta", &FirstPassageGreensFunction1DRad::seta )
	.def( "geta", &FirstPassageGreensFunction1DRad::geta )
	.def( "setr0", &FirstPassageGreensFunction1DRad::setr0 )
	.def( "getr0", &FirstPassageGreensFunction1DRad::getr0 )
	.def( "drawTime", &FirstPassageGreensFunction1DRad::drawTime )
	.def( "drawR", &FirstPassageGreensFunction1DRad::drawR )
	.def( "drawEventType", &FirstPassageGreensFunction1DRad::drawEventType )
	.def( "flux_tot", &FirstPassageGreensFunction1DRad::flux_tot )
	.def( "flux_rad", &FirstPassageGreensFunction1DRad::flux_rad )
	.def( "fluxRatioRadTot", &FirstPassageGreensFunction1DRad::fluxRatioRadTot )
	.def( "p_survival", &FirstPassageGreensFunction1DRad::p_survival )
	.def( "calcpcum", &FirstPassageGreensFunction1DRad::calcpcum )
	;
    
    class_<GreensFunction3DSym>("GreensFunction3DSym", init<Real>())
        .def( "getD", &GreensFunction3DSym::getD )
        .def( "drawTime", &GreensFunction3DSym::drawTime )
        .def( "drawR", &GreensFunction3DSym::drawR )
        .def( "p_r", &GreensFunction3DSym::p_r )
        .def( "ip_r", &GreensFunction3DSym::ip_r )
        .def( "dump", &GreensFunction3DSym::dump )
        ;

    class_<GreensFunction3DAbsSym>("GreensFunction3DAbsSym",
                                       init<Real, Real>())
        .def( "getD", &GreensFunction3DAbsSym::getD )
        .def( "geta", &GreensFunction3DAbsSym::geta )
        .def( "drawTime", &GreensFunction3DAbsSym::drawTime )
        .def( "drawR", &GreensFunction3DAbsSym::drawR )
        .def( "p_survival", &GreensFunction3DAbsSym::p_survival )
        .def( "p_int_r", &GreensFunction3DAbsSym::p_int_r )
        .def( "p_int_r_free", &GreensFunction3DAbsSym::p_int_r_free )
        //.def( "p_r_fourier", &GreensFunction3DAbsSym::p_r_fourier )
        ;

    class_<GreensFunction3DRadInf>("GreensFunction3DRadInf",
                                    init<Real, Real, Real, Real>())
        .def( "getD", &GreensFunction3DRadInf::getD )
        .def( "getkf", &GreensFunction3DRadInf::getkf )
        .def( "getSigma", &GreensFunction3DRadInf::getSigma )
        .def( "drawTime", &GreensFunction3DRadInf::drawTime )
        .def( "drawR", &GreensFunction3DRadInf::drawR )
        .def( "drawTheta", &GreensFunction3DRadInf::drawTheta )

//        .def( "p_tot", &GreensFunction3DRadInf::p_tot )
        .def( "p_free", &GreensFunction3DRadInf::p_free )
        .def( "ip_free", &GreensFunction3DRadInf::ip_free )
        .def( "p_corr", &GreensFunction3DRadInf::p_corr )
        .def( "ip_corr", &GreensFunction3DRadInf::ip_corr )
        .def( "p_survival", &GreensFunction3DRadInf::p_survival )
        .def( "p_int_r", &GreensFunction3DRadInf::p_int_r )
        .def( "p_theta", &GreensFunction3DRadInf::p_theta )
        .def( "ip_theta", &GreensFunction3DRadInf::ip_theta )

        .def( "dump", &GreensFunction3DRadInf::dump )
        ;

    class_<GreensFunction3D>("GreensFunction3D", init<Real, Real>())
        .def( "getD", &GreensFunction3D::getD )
        .def( "getkf", &GreensFunction3D::getkf )
        .def( "getSigma", &GreensFunction3D::getSigma )
        .def( "drawTime", &GreensFunction3D::drawTime )
        .def( "drawR", &GreensFunction3D::drawR )
        .def( "drawTheta", &GreensFunction3D::drawTheta )

        .def( "p_r", &GreensFunction3D::p_r )
        .def( "ip_r", &GreensFunction3D::ip_r )
        .def( "p_theta", &GreensFunction3D::p_theta )
        .def( "ip_theta", &GreensFunction3D::ip_theta )

        .def( "dump", &GreensFunction3D::dump )
        ;

    enum_<GreensFunction3DRadAbs::EventKind>("PairEventKind")
        .value( "IV_ESCAPE", GreensFunction3DRadAbs::IV_ESCAPE )
        .value( "IV_REACTION", GreensFunction3DRadAbs::IV_REACTION )
        ;

    class_<GreensFunction3DRadAbs>( "GreensFunction3DRadAbs",
                                            init<Real, Real, Real, Real, Real>() )
        .def( "geta", &GreensFunction3DRadAbs::geta )
        .def( "getD", &GreensFunction3DRadAbs::getD )
        .def( "getkf", &GreensFunction3DRadInf::getkf )
        .def( "getSigma", &GreensFunction3DRadInf::getSigma )
        .def( "drawTime", &GreensFunction3DRadAbs::drawTime )
        //.def( "drawTime2", &GreensFunction3DRadAbs::drawTime2 )
        .def( "drawEventType", &GreensFunction3DRadAbs::drawEventType )
        .def( "drawR", &GreensFunction3DRadAbs::drawR )
        .def( "drawTheta", &GreensFunction3DRadAbs::drawTheta )

        .def( "p_survival", &GreensFunction3DRadAbs::p_survival )
        .def( "dp_survival", &GreensFunction3DRadAbs::dp_survival )
        .def( "p_leaves", &GreensFunction3DRadAbs::p_leaves )
        .def( "p_leavea", &GreensFunction3DRadAbs::p_leavea )
        .def( "leaves", &GreensFunction3DRadAbs::leaves )
        .def( "leavea", &GreensFunction3DRadAbs::leavea )

        .def( "p_0", &GreensFunction3DRadAbs::p_0 )
        .def( "p_int_r", &GreensFunction3DRadAbs::p_int_r )
        .def( "p_int_r", &GreensFunction3DRadAbs::p_int_r )
        .def( "p_theta", &GreensFunction3DRadAbs::p_theta )
        .def( "ip_theta", &GreensFunction3DRadAbs::ip_theta )
        .def( "idp_theta", &GreensFunction3DRadAbs::idp_theta )

        .def( "f_alpha0", &GreensFunction3DRadAbs::f_alpha0 )
        .def( "alpha0_i", &GreensFunction3DRadAbs::alpha0_i )
        .def( "f_alpha", &GreensFunction3DRadAbs::f_alpha )
        .def( "f_alpha_aux", &GreensFunction3DRadAbs::f_alpha_aux )

        .def( "p_survival_i_exp", &GreensFunction3DRadAbs::p_survival_i_exp )
        .def( "p_survival_i_alpha", &GreensFunction3DRadAbs::p_survival_i_alpha )

        //.def( "guess_maxi", &GreensFunction3DRadAbs::guess_maxi )

        .def( "dump", &GreensFunction3DRadAbs::dump )

//        .def( "alpha_i", &GreensFunction3DRadAbs::alpha_i )
        ;

    class_<GreensFunction3DAbs>
        ( "GreensFunction3DAbs", init<Real, Real, Real>() ) 
        .def( "geta", &GreensFunction3DAbs::geta )
        .def( "getD", &GreensFunction3DAbs::getD )
        .def( "drawTime", &GreensFunction3DAbs::drawTime )
        .def( "drawR", &GreensFunction3DAbs::drawR )
        .def( "drawTheta", 
              &GreensFunction3DAbs::drawTheta )

        .def( "p_survival", 
              &GreensFunction3DAbs::p_survival )
        .def( "dp_survival", 
              &GreensFunction3DAbs::dp_survival )
        .def( "p_int_r", &GreensFunction3DAbs::p_int_r )
        .def( "p_theta", &GreensFunction3DAbs::p_theta )
        .def( "ip_theta", &GreensFunction3DAbs::ip_theta )
        .def( "idp_theta", 
              &GreensFunction3DAbs::idp_theta )

        .def( "dump", &GreensFunction3DAbs::dump )
        ;
}
