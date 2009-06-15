#if !defined( __PAIRGREENSFUNCTION_HPP )
#define __PAIRGREENSFUNCTION_HPP

#include <vector>
#include <boost/multi_array.hpp>

#include "Defs.hpp"


enum EventType
{
    REACTION = 0,
    ESCAPE = 1
};



class PairGreensFunction
{

public:

  PairGreensFunction( const Real D, const Real kf, const Real Sigma )
    :
    D( D ),
    kf( kf ),
    Sigma( Sigma )
  {
    ;
  }

  ~PairGreensFunction()
  {
    ;
  }

  const Real getD() const
  {
    return this->D;
  }

  const Real getkf() const
  {
    return this->kf;
  }

  const Real getSigma() const
  {
    return this->Sigma;
  }



/*
  const Real drawTime( const Real rnd, const Real r0,
                               const Real maxt ) const = 0;

  const Real drawR( const Real rnd, 
                            const Real r0, 
                            const Real t ) const = 0;

  const Real drawTheta( const Real rnd,
                                const Real r, 
                                const Real r0, 
                                const Real t ) const = 0;
*/

private:

  const Real D;
  const Real kf;
  const Real Sigma;

};



#endif // __PAIRGREENSFUNCTION_HPP
