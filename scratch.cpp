

///
/// doesn't run fast, but if bessel_sequence uses adaptive runge kutta
/// it should show excellent performance.
//
const Real 
FirstPassagePairGreensFunction::p_n( const Integer n,
				     const Real r,
				     const Real r0, 
				     const Real t ) const
{
    Real p( 0.0 );

    const Real a( geta() );
    const Real sigma( getSigma() );

    Real realn( 0.5 + n );

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    RealVector::size_type mini( 0 );
    RealVector::size_type maxi( 4 );

    RealVector alphaTable;
    RealVector jas1Table;
    RealVector jaaTable;
    RealVector jarTable;
    RealVector jar0Table;

    RealVector pTable;

    const unsigned int offset( alphaOffset( n ) );    

    while( true )
    {
        for( RealVector::size_type i( 0 ); i < maxi-mini; ++i )
        {
            alphaTable.push_back( alpha_i( i+mini+offset, n, solver ) );
            jas1Table.push_back( alphaTable[i] * sigma );
            jaaTable.push_back( alphaTable[i] * a );
            jarTable.push_back( alphaTable[i] * r );
            jar0Table.push_back( alphaTable[i] * r0 );
        }
        
        RealVector jas2Table( jas1Table );

        gsl_sf_bessel_sequence_Jnu_e( realn, GSL_PREC_SINGLE,
                                      jas1Table.size(), &jas1Table[0] );
        gsl_sf_bessel_sequence_Jnu_e( realn+1.0, GSL_PREC_SINGLE,
                                      jas2Table.size(), &jas2Table[0] );
        gsl_sf_bessel_sequence_Jnu_e( realn, GSL_PREC_SINGLE,
                                      jaaTable.size(), &jaaTable[0] );
        gsl_sf_bessel_sequence_Jnu_e( realn, GSL_PREC_SINGLE,
                                      jarTable.size(), &jarTable[0] );
        gsl_sf_bessel_sequence_Jnu_e( realn, GSL_PREC_SINGLE,
                                      jar0Table.size(), &jar0Table[0] );

        
        for( RealVector::size_type i( 0 ); i < maxi-mini; ++i )
        {
            const Real p_i( p_n_alpha( alphaTable[i], n, r, r0, t, jas1Table[i], 
                                       jas2Table[i], jaaTable[i], jarTable[i],
                                       jar0Table[i] ) );
            
            p += p_i;
            pTable.push_back( p_i );
        }

        if( fabs( p * TOLERANCE ) > fabs( pTable.back() ) || maxi >= MAX_ALPHA_SEQ )
        {
//            printf("%d %g %g\n",maxi, p, pTable.back() );
            break;
        }

        alphaTable.clear();
        jas1Table.clear();
        jas2Table.clear();
        jaaTable.clear();
        jarTable.clear();
        jar0Table.clear();

        mini = maxi;
        maxi *= 2;
    }

    //p = std::accumulate( pTable.begin(), pTable.end(), 0.0 );
    return p;
//    return pTable.back();
}


#if 0
const Real PlainPairGreensFunction::drawTheta( const Real rnd,
					       const Real r, 
					       const Real r0, 
					       const Real t )
{
  const Real thetaStep( getThetaStep( r, r0, t ) );
  
  const Index tableSize( getTableSizeOfTheta() );

  RealVector pTable( tableSize );

  Real p_prev( 0.0 ); // this->p_tot( r, r0, 0, t ) );
  for( Index thetaIndex( 1 ); thetaIndex < tableSize;
       ++thetaIndex )
    {
      const Real theta( thetaStep * thetaIndex );

      const Real p( this->p_tot( r, r0, theta, t ) );

      pTable[thetaIndex] = pTable[thetaIndex-1] + 
	( p_prev + p ) * 0.5 * thetaStep;

      p_prev = p;
    }

  //  const Real p_irr_r( this->p_irr_radial( r, t, r0 ) );
  //  printf( "%g %g\n", pTable[pTable.size()-1], p_irr_r );



  return 0.0;
}
#endif// 0


#if 0
const Real PlainPairGreensFunction::drawR( const Real rnd, 
					   const Real r0, 
					   const Real t )
{
  const Real tScale( getScaleOfT( t ) );
  const Index tIndex( static_cast<Index>
		      ( round( tScale * getTableSizeOfT() ) ) );
  const Real tOnTable( getTByScale( (Real)tIndex / getTableSizeOfT() ) );

  const Real r0Scale( getScaleOfR0( r0, tOnTable ) );
  const Index r0Index( static_cast<Index>
		       ( round( r0Scale * getTableSizeOfR0() ) ) );
  const Real r0OnTable( getR0ByScale( (Real)r0Index / getTableSizeOfR0(), 
				      t ) );
  RealVector ans( getTableSizeOfR() );

  // ans[0] = 0.0;  // no need to do this; std::vector is filled with zero.
  for( Index rIndex( 1 ); rIndex < getTableSizeOfR(); ++rIndex )
    {

      /*

      const Real thetaStep( getThetaStep( 0,0,0 ) );

      Real p_cum_theta( 0.0 );
      Real value_theta_prev( lookupP_totTable( rIndex, r0Index, 0, tIndex )  );
      for( Index thetaIndex( 1 ); thetaIndex <= getTableSizeOfTheta(); 
	   ++thetaIndex )
	{
	  const Real value( lookupP_totTable( rIndex, r0Index, 
					      thetaIndex, tIndex ) );
	  p_cum_theta += ( value_theta_prev + value ) * 0.5 * thetaStep;

	  value_theta_prev = value;
	}
      */

      const Real p( p_irr_radial( getRByScale( (Real)rIndex / 
					       getTableSizeOfR(), t ),
				  tOnTable, r0OnTable ) );
      ans[rIndex] = ans[rIndex-1] + p;
    }

  const Real targetPoint( rnd * *(ans.end()) );

  const size_t i( gsl_interp_bsearch( ans.data(), targetPoint, 
				      0, ans.size()-1 ) );

  const Real lowR( getRByScale( (Real)i / getTableSizeOfR(), t ) );
  const Real highR( getRByScale( (Real)(i+1) / getTableSizeOfR(), t ) );

  const Real r( lowR * 
		( 1 + ( targetPoint - ans[i] ) / ( ans[i+1] - ans[i] ) ) );

  //printf("%g %g %g %g\n",ans[i],ans[i+1],targetPoint,r);

  //  printf("%g %g\n", r0, getR0ByScale( r0Scale, t ) );

  return r;
}
#endif //  0



  const Index getTableSizeOfT() const
  {
    return 20;
  }

  const Index getTableSizeOfR() const
  {
    return 40;
  }

  const Index getTableSizeOfR0() const
  {
    return 40;
  }

  const Index getTableSizeOfTheta() const
  {
    return 1000;
  }

  const Real getMinT() const
  {
    return sqrt( getD() );
  }

  const Real getMaxT() const
  {
    return getMinT() * 1e3;
  }

  const Real getScaleOfT( const double t ) const
  {
    const Real tmin( getMinT() );
    const Real tmax( getMaxT() );  
    
    const Real scale( log( t / tmin ) / log( tmax/tmin ) );

    return scale;
  }

  const Real getTByScale( const double scale ) const
  {
    const Real tmin( getMinT() );
    const Real tmax( getMaxT() );
    
    const Real tstep( log( tmax / tmin ) ); 

    return tmin * exp( tstep * scale );
  }

  const Real getMinR0() const
  {
    return getSigma();
  }

  const Real getMaxR0( const Real t ) const
  {
    // [ Sigma, Sigma + H sqrt( 6 D t ) ]
    return this->H * sqrt( 6.0 * getD() * t );
  }

  const Real getScaleOfR0( const Real r0, const Real t ) const
  {
    const Real r0min( getMinR0() );
    const Real r0max( getMaxR0( t ) );  
    const Real r0range( r0max - r0min );

    return ( r0 - r0min ) / r0range;
  }

  const Real getR0ByScale( const Real scale, const Real t ) const
  {
    const Real r0min( getMinR0() );
    const Real r0max( getMaxR0( t ) );  
    const Real r0range( r0max - r0min );

    return r0min + r0range * scale;
  }

  const Real getRByScale( const Real scale, const Real t ) const
  {
    return getMinR() + ( getMaxR( t ) - getMinR() ) * scale;
  }

  const Real getMinTheta() const
  {
    return 0.001;
    //return 0.0;
  }

  const Real getMaxTheta( const Real r, const Real r0, const Real t ) const
  {
    return M_PI; //* 0.5;

    /*
    double value( 6.0 * pow( (r0+r)/sqrt(6.0 * getD() * t), -0.4 ) - 1.0 );
    if( value > M_PI )
      {
	value = M_PI;
      }

    return value;
    */
  }

  const Real getThetaStep( const Real r, const Real r0, const Real t ) const
  {
    return ( getMaxTheta( r, r0, t ) - getMinTheta() ) / getTableSizeOfTheta();
  }

  const Real getThetaByScale( const Real scale, const Real r,
			      const Real r0, const Real t ) const
  {
    return getMinTheta() + 
      ( getMaxTheta( r, r0, t ) - getMinTheta() ) * scale;
  }


#if 0
const Real 
PlainPairGreensFunction::lookupP_totTable( const Index rIndex, 
					   const Index r0Index,
					   const Index thetaIndex, 
					   const Index tIndex )
{
  Real value( p_totTable[rIndex][r0Index][thetaIndex][tIndex] );

  if( value < 0.0 ) // not in cache
    {
      const Real t( getTByScale( Real( tIndex ) / getTableSizeOfT() ) );
      const Real r( getRByScale( Real( rIndex ) / getTableSizeOfR(), t ) );
      const Real r0( getR0ByScale( Real( r0Index ) / getTableSizeOfR0(), t ) );
      const Real theta( getThetaByScale( Real( thetaIndex ) / 
					 getTableSizeOfTheta(), t, r0, r ) );
      
      //  printf("table %g %g %g %g\n",r, r0, theta,t);
      
      value = p_tot( r, r0, theta, t );
      
      p_totTable[rIndex][r0Index][thetaIndex][tIndex] = value;
      
    }

  return value;
}
#endif // 0



  HalfOrderBesselGenerator us( u * params->Sigma, order, order );
  HalfOrderBesselGenerator ur( u * params->r, order, order );
  HalfOrderBesselGenerator ur0( u * params->r0, order, order );
      const Real js( us.j( order ) );
      const Real ys( us.y( order ) );
      const Real jps( us.jp( order ) );
      const Real yps( us.yp( order ) );
      const Real jr( ur.j( order ) );
      const Real yr( ur.y( order ) );
      const Real jr0( ur0.j( order ) );
      const Real yr0( ur0.y( order ) );



      //Real jr,yr,jr0,yr0,js,ys,jps,yps,tmp_;
      //    bessjy(u*r,order+0.5,&jr,&yr,&tmp_,&tmp_);
      //      bessjy(u*r0,order+0.5,&jr0,&yr0,&tmp_,&tmp_);
      //      bessjy(u*Sigma,order+0.5,&js,&ys,&jps,&yps);






const Real 
PlainPairGreensFunction::intt_p_irr_radial( const Real r, 
					    const Real r0, 
					    const Real t ) const
{
  static const Real M_1_SQRTPI( M_2_SQRTPI * 0.5 );

  const Real rmr0( r - r0 );
  const Real rpr0m2Sigma( r + r0 - 2.0 * getSigma() );
  const Real Dt4( 4.0 * getD() * t );
  const Real sqrtDt4( sqrt( Dt4 ) );
  const Real sqrtD( sqrt( getD() ) );

  const Real alpha( this->alpha );

  const Real num1( M_1_SQRTPI * sqrtD *
		   ( exp( - rmr0 * rmr0 / Dt4 ) -
		     exp( - rpr0m2Sigma * rpr0m2Sigma / Dt4 ) ) );
  const Real num2( sqrtD * 0.5 *
		   rmr0 * erf( rmr0 / sqrtDt4 ) );
  const Real num3( ( sqrtD * 0.5 / alpha ) *
		   ( alpha * rpr0m2Sigma + 2.0 * sqrtD ) *
		   erf( rpr0m2Sigma / sqrtDt4 ) );

  const Real num4( ( 1.0 / alpha ) * 
  		   W( rpr0m2Sigma / sqrtDt4, alpha * sqrt( t ) ) );

  const Real result( ( num1 + num2 - num3 - num4 ) / 
		     ( 4.0 * sqrtD * M_PI * r * r0 ) );

  //  printf("n %g %g %g %g %g %g\n", num1,num2,num3,num4,( 4.0 * sqrtD * M_PI * r * r0 ),result );


  return result; //  * 4.0 * M_PI * r * r;
}





const Real PlainPairGreensFunction::drawTime( const Real rnd, 
					      const Real r0, 
					      const Real maxt ) const
{
  assert( rnd <= 1.0 && rnd >= 0.0 );

  {
    const Real maxp( p_survival( maxt, r0 ) );

    if( rnd >= maxp )
      {
	return INFINITY;
      }
  }

  p_survival_params params = { this, r0, rnd };
  gsl_function F = 
    {
      reinterpret_cast<typeof(F.function)>( &p_survival_F ),
      &params 
    };

  // This initial guess must be smaller than the answer, or
  // newton-type iteration can undershoot into tsqrt < 0, where
  // value of p_survival is undefined.
  const Real initialGuess( 1e-100 );

  //  const gsl_root_fdfsolver_type* solverType( gsl_root_fdfsolver_steffenson );
  const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
  gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
  gsl_root_fsolver_set( solver, &F, initialGuess, maxt );

  const Index maxIter( 100 );

  Real tsqrt( initialGuess );
  Index i( 0 );
  while( true )
    {
      gsl_root_fsolver_iterate( solver );

      tsqrt = gsl_root_fsolver_root( solver );
      const Real x_lo( gsl_root_fsolver_x_lower( solver ) );
      const Real x_hi( gsl_root_fsolver_x_upper( solver ) );

      //      int status( gsl_root_test_delta( tsqrt, tsqrt_prev, 0.0, 
      //				       1e-12 ) );
      int status = gsl_root_test_interval( x_lo, x_hi, 1e-30, 1e-18 );
      if( status == GSL_CONTINUE )
	{
	  if( i >= maxIter )
	    {
	      gsl_root_fsolver_free( solver );
	      std::cerr << "drawTime: failed to converge" << std::endl;
	      throw std::exception();
	    }
	}
      else
	{
	  break;
	}

      ++i;
    }
  
  gsl_root_fsolver_free( solver );

  return gsl_pow_2( tsqrt );
} 



const Real PlainPairGreensFunction::drawTime( const Real rnd, 
					      const Real r0, 
					      const Real maxt ) const
{
  assert( rnd <= 1.0 && rnd >= 0.0 );

  {
    const Real maxp( p_survival( maxt, r0 ) );

    if( rnd >= maxp )
      {
	return INFINITY;
      }
  }

  p_survival_params params = { this, r0, rnd };
  gsl_function_fdf FDF = 
    {
      reinterpret_cast<typeof(FDF.f)>( &p_survival_F ),
      reinterpret_cast<typeof(FDF.df)>( &p_survival_deriv_F ),
      reinterpret_cast<typeof(FDF.fdf)>( &p_survival_fdf_F ),
      &params 
    };

  // This initial guess must be smaller than the answer, or
  // newton-type iteration can undershoot into tsqrt < 0, where
  // value of p_survival is undefined.
  const Real initialGuess( 1e-100 );

  //  const gsl_root_fdfsolver_type* solverType( gsl_root_fdfsolver_steffenson );
  const gsl_root_fdfsolver_type* solverType( gsl_root_fdfsolver_secant );
  gsl_root_fdfsolver* solver( gsl_root_fdfsolver_alloc( solverType ) );
  gsl_root_fdfsolver_set( solver, &FDF, initialGuess );

  const Index maxIter( 100 );

  Real tsqrt( initialGuess );
  Index i( 0 );
  while( true )
    {
      gsl_root_fdfsolver_iterate( solver );
      const Real tsqrt_prev( tsqrt );
      tsqrt = gsl_root_fdfsolver_root( solver );

      int status( gsl_root_test_delta( tsqrt, tsqrt_prev, 0.0, 
				       1e-12 ) );

      if( status == GSL_CONTINUE )
	{
	  if( i >= maxIter )
	    {
	      gsl_root_fdfsolver_free( solver );
	      std::cerr << "drawTime: failed to converge" << std::endl;
	      throw std::exception();
	    }
	}
      else
	{
	  break;
	}

      ++i;
    }
  
  gsl_root_fdfsolver_free( solver );

  return gsl_pow_2( tsqrt );
} 



const Real PlainPairGreensFunction::drawTime( const Real rnd, 
					      const Real r0, 
					      const Real maxt ) const
{
  assert( rnd <= 1.0 && rnd >= 0.0 );

  gsl_function_fdf FDF;
  p_survival_params params = {this, r0, rnd};
     
      
  FDF.f = reinterpret_cast<typeof(FDF.f)>( &p_survival_F );
  FDF.df = reinterpret_cast<typeof(FDF.df)>( &p_survival_deriv_F );
  FDF.fdf = reinterpret_cast<typeof(FDF.fdf)>( &p_survival_fdf_F );
  FDF.params = &params;
  

  Real x( 7e-100 );

  const gsl_root_fdfsolver_type* T( gsl_root_fdfsolver_newton );
  gsl_root_fdfsolver* s( gsl_root_fdfsolver_alloc(T) );
  gsl_root_fdfsolver_set (s, &FDF, x);
  printf ("%-5s %10s %10s %10s\n",
	  "iter", "root", "err", "err(est)");
  
  int status;
  

  int iter = 0, max_iter = 100;
  do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      Real x0( x );
      x = gsl_root_fdfsolver_root (s);

      status = gsl_root_test_delta (x, x0, 0, 1e-16);
      
      if (status == GSL_SUCCESS)
	printf ("Converged:\n");
      
      printf ("%5d %10.7f %10.7f\n",
	      iter, x, x - x0);
    }
  while (status == GSL_CONTINUE && iter < max_iter);
  
  gsl_root_fdfsolver_free( s );
  return x*x;
} 


const Real PlainPairGreensFunction::drawTime( const Real rnd, 
					      const Real r0, 
					      const Real maxt ) const
{
  assert( rnd <= 1.0 && rnd >= 0.0 );

  gsl_function F;
  p_survival_params params = {this, r0, rnd};
     
  F.function = reinterpret_cast<typeof(F.function)>( &p_survival_F );
  F.params = &params;

  const gsl_root_fsolver_type* T( gsl_root_fsolver_brent );
  gsl_root_fsolver* s( gsl_root_fsolver_alloc(T) );
  gsl_root_fsolver_set( s, &F, 1e-100, maxt );
     
  int status;

  int iter = 0, max_iter = 100;
  Real t;
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      t = gsl_root_fsolver_root(s);
      const Real x_lo( gsl_root_fsolver_x_lower(s) );
      const Real x_hi( gsl_root_fsolver_x_upper(s) );
      status = gsl_root_test_interval( x_lo, x_hi, 0, 1e-16 );
      
//      printf ("%5d %10.7f %10.7f\n",   iter, x_lo, x_hi);
    }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free( s );

  return t*t;
} 



/*
const Real 
PlainPairGreensFunction::p_survival_deriv( const Real t, 
					   const Real r0 ) const
{
  Real deriv;

  const Real Sigma( this->getSigma() );
  const Real D( this->getD() );
  const Real alpha( this->getalpha() );
  const Real kD( this->getkD() );
  const Real kf( this->getkf() );

  const Real r0_m_Sigma_over_sqrt4Dt( ( r0 - Sigma ) / sqrt( 4.0 * D * t ) );

  const Real sqrtPit( sqrt( M_PI * t ) );

  const Real num1( exp( - gsl_pow_2( r0_m_Sigma_over_sqrt4Dt ) ) );
  const Real num2( alpha * sqrtPit * 
  W( r0_m_Sigma_over_sqrt4Dt, alpha * sqrt( t ) ) );



  const Real factor( alpha * kf * Sigma / ( sqrtPit * r0 * ( kf + kD ) ) );

  deriv = ( num2 - num1 ) * factor;

  return deriv;
}

void
PlainPairGreensFunction::p_survival_fdf( const Real t, 
					 const Real r0,
					 Real* const f, Real* const df ) const
{
  const Real kD( this->getkD() );
  const Real kf( this->getkf() );
  const Real Sigma( this->getSigma() );
  const Real D( this->getD() );
  const Real alpha( this->getalpha() );
  const Real sqrtPit( sqrt( M_PI * t ) );

  const Real r0_m_Sigma_over_sqrt4Dt( ( r0 - Sigma ) / sqrt( 4.0 * D * t ) );
  const Real W( W( r0_m_Sigma_over_sqrt4Dt, alpha * sqrt( t ) ) );
  
  const Real Sigmakf_over_r0_kf_p_kD( Sigma * kf / ( r0 * ( kf + kD ) ) );

  *f = 1.0 - Sigmakf_over_r0_kf_p_kD * ( erfc( r0_m_Sigma_over_sqrt4Dt ) - W );

  const Real dfnum1( exp( - gsl_pow_2( r0_m_Sigma_over_sqrt4Dt ) ) );
  const Real dfnum2( alpha * sqrtPit * W );
  const Real dffactor( ( alpha / sqrtPit ) * Sigmakf_over_r0_kf_p_kD );

  *df = ( dfnum2 - dfnum1 ) * dffactor;
}
*/




const Real 
FirstPassagePairGreensFunction::asratio( const Real alpha,
					 const Real r0 ) const
{
    const Real a( geta() );
    const Real D( getD() );
    const Real sigma( getSigma() );
    const Real h( geth() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real sigmasq( sigma * sigma );
    const Real alphasq( alpha * alpha );

    const Real angle_a( alpha * ( a - sigma ) );
    Real sin_a;
    Real cos_a;
    sincos( angle_a, &sin_a, &cos_a );
    const Real num( - a * ( ( hsigma_p_1 ) * cos_a -
			    sigma * alpha * sin_a ) );
		      

    const Real den( h * sigmasq );


    const Real result( num / den );

    return result;
}

const Real
FirstPassagePairGreensFunction::
f_alpha0_aux_df_F( const Real alpha,
		   const f_alpha0_aux_params* const params )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 

    return gf->f_alpha0_aux_df( alpha );
}


void
FirstPassagePairGreensFunction::
f_alpha0_aux_fdf_F( const Real alpha,
		    const f_alpha0_aux_params* const params,
		    Real* const f, Real* const df )
{
    const FirstPassagePairGreensFunction* const gf( params->gf ); 
    const Real value( params->value );

    *f = gf->f_alpha0_aux( alpha ) - value;
    *df = gf->f_alpha0_aux_df( alpha );
}

const Real 
FirstPassagePairGreensFunction::
f_alpha0_aux_df( const Real alpha ) const

{
    const Real a( geta() );
    const Real sigma( getSigma() );
    const Real hsigma_p_1( this->hsigma_p_1 );

    const Real alphasq( alpha * alpha );

    const Real term1( alpha - sigma );
    const Real term2( hsigma_p_1 / 
		      ( sigma * alphasq * 
			( 1.0 + ( hsigma_p_1 * hsigma_p_1 / 
				  ( sigma * sigma * alphasq ) ) ) ) );

    const Real result( term1 + term2 );

    return result;
}
