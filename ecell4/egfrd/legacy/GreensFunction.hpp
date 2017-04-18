#if !defined( EGFRD_GREENSFUNCTION_HPP )
#define EGFRD_GREENSFUNCTION_HPP

#include "Defs.hpp"

class GreensFunction
{
public:
    enum EventKind
    {
        IV_ESCAPE,
        IV_REACTION
    };

public:
    GreensFunction( const Real D )
      : D( D ) {}
  
    ~GreensFunction() {}
  
    Real getD() const
    {
        return this->D;
    }

protected:
    const Real D;
};

#endif // EGFRD_GREENSFUNCTION_HPP
