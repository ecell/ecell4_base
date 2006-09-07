#if !defined( __PAIRGREENSFUNCTION_HPP )
#define __PAIRGREENSFUNCTION_HPP

#include <vector>
#include <boost/multi_array.hpp>


typedef double Real;
typedef int Int;
typedef size_t Index;


typedef std::vector< Real > RealVector;
//typedef boost::multi_array< Real, 1, boost::pool_allocator<Real> > 
//RealArray;
typedef boost::multi_array<Real, 2>
Real2DArray;
typedef boost::multi_array<Real, 3>
Real3DArray;
typedef boost::multi_array<Real, 4>
Real4DArray;


class PairGreensFunction
{

public:

  PairGreensFunction( const Real D, const Real kf, const Real Sigma )
    :
    D( D ),
    kf( kf ),
    Sigma( Sigma )
  {
    ;
  }

  virtual ~PairGreensFunction()
  {
    ;
  }

  const Real getD() const
  {
    return this->D;
  }

  const Real getkf() const
  {
    return this->kf;
  }

  const Real getSigma() const
  {
    return this->Sigma;
  }

  virtual const Real p_tot( const Real r, const Real r0, 
			    const Real theta, const Real time ) const = 0;

  //  virtual const Real p_survival( const Real t, const Real r0 ) = 0;

  virtual const Real drawNextReactionTime( const Real rnd, const Real r0 ) = 0;

  virtual const Real drawR( const Real r, 
			    const Real r0, 
			    const Real t ) const = 0;

  virtual const Real drawTheta( const Real theta, 
				const Real r, 
				const Real r0, 
				const Real t ) const = 0;


private:

  const Real D;
  const Real kf;
  const Real Sigma;

};



#endif // __PAIRGREENSFUNCTION_HPP
