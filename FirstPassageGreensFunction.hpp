#if !defined( __FIRSTPASSAGEGREENSFUNCTION_HPP )
#define __FIRSTPASSAGEGREENSFUNCTION_HPP

#include <vector>
#include <boost/multi_array.hpp>

#include "Defs.hpp"


class FirstPassageGreensFunction
{

public:

    FirstPassageGreensFunction( const Real D )
	:
	D( D ),
        a( 0.0 )
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

    void seta( const Real a )
    {
        THROW_UNLESS( std::invalid_argument, a >= 0.0 );
        this->a = a;
    }

    const Real geta() const
    {
        return this->a;
    }

    const Real p_survival( const Real t ) const; 

    const Real drawTime( const Real rnd ) const;

    const Real drawR( const Real rnd, const Real t ) const;

    const Real p_r_int( const Real r, const Real t ) const;
    const Real p_free_int( const Real r, const Real t ) const;

    const Real p_r_fourier( const Real r, const Real t ) const;

private:

    struct p_survival_params
    {
	const FirstPassageGreensFunction* const gf;
	const Real rnd;
    };

    static const Real p_survival_F( const Real t, 
				    const p_survival_params* params );

    struct p_r_params
    {
	const FirstPassageGreensFunction* const gf;
	const Real t;
	const Real St;
	const Real rnd;
    };

    static const Real p_r_F( const Real r, 
			     const p_r_params* params );


private:

    static const Real CUTOFF = 1e-10;

    const Real D;

    Real a;
};



#endif // __PAIRGREENSFUNCTION_HPP
