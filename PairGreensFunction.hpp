#if !defined( __PAIRGREENSFUNCTION_HPP )
#define __PAIRGREENSFUNCTION_HPP

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

class PairGreensFunction
{
public:
    PairGreensFunction(Real D, Real kf, Real Sigma)
      : D(D), kf(kf), Sigma(Sigma) {}
    
    ~PairGreensFunction() {}
    
    Real getD() const
    {
        return this->D;
    }
    
    Real getkf() const
    {
        return this->kf;
    }
    
    Real getSigma() const
    {
        return this->Sigma;
    }

private:
  const Real D;
  const Real kf;
  const Real Sigma;
};

#endif /* __PAIRGREENSFUNCTION_HPP */
