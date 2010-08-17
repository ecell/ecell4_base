#if !defined( __FIRSTPASSAGEGREENSFUNCTION_HPP)
#define __FIRSTPASSAGEGREENSFUNCTION_HPP

#include "Defs.hpp"

#include "Logger.hpp"

class GreensFunction3DAbsSym
{
public:
    GreensFunction3DAbsSym(Real D, Real a)
        : D( D), a( a) {}

    ~GreensFunction3DAbsSym() {}

    Real getD() const
    {
        return this->D;
    }

    Real geta() const
    {
        return this->a;
    }

    Real p_survival(Real t) const; 

    Real drawTime(Real rnd) const;

    Real drawR(Real rnd, Real t) const;

    Real p_int_r(Real r, Real t) const;
    Real p_int_r_free(Real r, Real t) const;

    Real p_r_fourier(Real r, Real t) const;

    std::string dump() const;

private:
    static Real ellipticTheta4Zero(Real q);

private:

    static const Real CUTOFF = 1e-10;

    // H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    // 5.6: ~1e-8, 6.0: ~1e-9
    static const Real CUTOFF_H = 6.0;

    const Real D;
    const Real a;

    static Logger& log_;
};



#endif // __PAIRGREENSFUNCTION_HPP
