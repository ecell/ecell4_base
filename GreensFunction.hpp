#if !defined( __GREENSFUNCTION_HPP )
#define __GREENSFUNCTION_HPP

#include "Defs.hpp"

enum EventType
{
    SINGLE_REACTION = 10,
    SINGLE_ESCAPE = 11,

    COM_ESCAPE = 11, // Same as SINGLE_ESCAPE

    IV_EVENT = 12,
    IV_ESCAPE = 13,
    IV_REACTION = 14,

    IV_INTERACTION = 15,

    BURST = 16,

    MULTI_ESCAPE = 17,
    MULTI_REACTION = 18
};

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
