#if !defined( __GREENSFUNCTION_HPP )
#define __GREENSFUNCTION_HPP

#include "Defs.hpp"

class GreensFunction
{

public:

  GreensFunction( const Real D )
    :
    D( D )
  {
    ;
  }

  ~GreensFunction()
  {
    ;
  }

  const Real getD() const
  {
    return this->D;
  }

private:

  const Real D;

};



#endif // __GREENSFUNCTION_HPP
