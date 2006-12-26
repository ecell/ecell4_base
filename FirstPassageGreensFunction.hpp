#if !defined( __FIRSTPASSAGEGREENSFUNCTION_HPP )
#define __FIRSTPASSAGEGREENSFUNCTION_HPP

#include <vector>
#include <boost/multi_array.hpp>

#include "Defs.hpp"


class FirstPassageGreensFunction
{

public:

  FirstPassageGreensFunction( const Real D, const Real a )
    :
    D( D ),
    a( a )
  {
    ;
  }

  virtual ~FirstPassageGreensFunction()
  {
    ;
  }

  const Real getD() const
  {
    return this->D;
  }

  const Real geta() const
  {
    return this->a;
  }


  const Real p_survival( const Real time ) const; 

  const Real drawExitTime( const Real rnd, const Real r,
				   const Real time ) const;

  const Real drawR( const Real rnd, const Real r, const Real t ) const;

  const Real p_r_int( const Real r, const Real time ) const;
  const Real p_free_int( const Real r, const Real time ) const;

  const Real p_r_fourier( const Real r, const Real time ) const;




private:

  const Real D;
  const Real a;
};



#endif // __PAIRGREENSFUNCTION_HPP
