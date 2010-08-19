#if !defined( __PAIRGREENSFUNCTION_HPP )
#define __PAIRGREENSFUNCTION_HPP

#include "Defs.hpp"
#include "GreensFunction.hpp"

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
