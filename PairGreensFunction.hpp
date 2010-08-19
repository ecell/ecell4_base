#if !defined( __PAIRGREENSFUNCTION_HPP )
#define __PAIRGREENSFUNCTION_HPP

#include "Defs.hpp"
#include "GreensFunction.hpp"

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

class PairGreensFunction: public GreensFunction
{
public:
    PairGreensFunction(Real D, Real kf, Real r0, Real Sigma)
      : GreensFunction(D), kf(kf), r0(r0), Sigma(Sigma) {}
    
    ~PairGreensFunction() {}
    
    Real getkf() const
    {
        return this->kf;
    }
    
    Real getSigma() const
    {
        return this->Sigma;
    }

    Real getr0() const
    {
        return this->r0;
    }

protected:
  const Real kf;
  const Real r0;
  const Real Sigma;
};

#endif /* __PAIRGREENSFUNCTION_HPP */
