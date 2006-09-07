

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
