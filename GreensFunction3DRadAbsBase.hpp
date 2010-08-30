#ifndef FIRST_PASSAGE_PAIR_GREENS_FUNCTION_BASE_HPP
#define FIRST_PASSAGE_PAIR_GREENS_FUNCTION_BASE_HPP

#include "PairGreensFunction.hpp"

class GreensFunction3DRadAbsBase: public PairGreensFunction
{
public:
    enum EventKind
    {
        IV_ESCAPE,
        IV_REACTION
    };

public:
    GreensFunction3DRadAbsBase(Real D, Real kf, Real r0, Real Sigma)
        : PairGreensFunction(D, kf, r0, Sigma) {}

    virtual ~GreensFunction3DRadAbsBase();

    virtual Real drawTime(Real rnd) const = 0;

    virtual EventKind drawEventType(Real rnd, Real t) const = 0;
    
    virtual Real drawR(Real rnd, Real t) const = 0;
    
    virtual Real drawTheta(Real rnd, Real r, Real t) const = 0;
}; 

#endif /* FIRST_PASSAGE_PAIR_GREENS_FUNCTION_BASE_HPP */
