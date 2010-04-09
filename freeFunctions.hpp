#ifndef FREE_FUNTIONS_HPP
#define FREE_FUNTIONS_HPP

#include "Defs.hpp"

const Real expxsq_erfc( const Real x );
const Real W( const Real a, const Real b );

const Real p_irr( const Real r, const Real r0, const Real t,
                  const Real kf, const Real D, const Real sigma );

const Real __p_irr( const Real r, const Real r0, const Real t,
                    const Real kf, const Real D, const Real sigma, 
                    const Real alpha );

const Real p_free( const Real r, const Real r0, 
                   const Real theta, const Real t );

const Real p_survival_irr( const Real t, const Real r0,
                  const Real kf, const Real D, const Real sigma );

const Real __p_reaction_irr( const Real t, const Real r0,
                             const Real kf, const Real D, const Real sigma,
                             const Real alpha, const Real kD );

const Real __p_reaction_irr_t_inf( const Real r0, const Real kf, 
                                   const Real sigma, const Real kD );

const Real p_survival_nocollision( const Real t,
                                   const Real r0,
                                   const Real D, const Real a );

const Real dp_survival_nocollision( const Real t,
                                    const Real r0,
                                    const Real D, const Real a );

const Real p_theta_free( const Real theta, const Real r, const Real r0, 
                         const Real t, const Real D );

const Real ip_theta_free( const Real theta, const Real r, const Real r0,
                          const Real t, const Real D );

const Real g_bd( const Real r0, const Real sigma, const Real t, const Real D );

const Real I_bd( const Real sigma, const Real t, const Real D );

const Real I_bd_r( const Real r, const Real sigma, const Real t, const Real D );

const Real drawR_gbd( const Real rnd, const Real sigma, 
                      const Real t, const Real D );

#endif /* FREE_FUNTIONS_HPP */
