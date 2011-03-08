#if !defined( __PAIRGREENSFUNCTION_HPP )
#define __PAIRGREENSFUNCTION_HPP

#include "Defs.hpp"
#include <string>
#include "GreensFunction.hpp"

class PairGreensFunction: public GreensFunction
{
public:
    PairGreensFunction(Real D, Real kf, Real r0, Real Sigma)
      : GreensFunction(D), kf(kf), r0(r0), Sigma(Sigma) {}
    
    virtual ~PairGreensFunction() {}
    
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

    Real getr0() const
    {
        return this->r0;
    }

    virtual std::string dump() const = 0;

    virtual const char* getName() const = 0;

    virtual Real drawTime(Real rnd) const = 0;
    
    virtual Real drawR(Real rnd, Real t) const = 0;
    
    virtual Real drawTheta(Real rnd, Real r, Real t) const = 0;

protected:
  const Real kf;
  const Real r0;
  const Real Sigma;
};

#endif /* __PAIRGREENSFUNCTION_HPP */
