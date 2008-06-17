#include "SphericalBesselTable.hpp"

#include "SphericalBesselGenerator.hpp"


const SphericalBesselGenerator& getSphericalBesselGenerator()
{
    static const SphericalBesselGenerator sphericalBesselGenerator;
    return sphericalBesselGenerator;
}

static const Real interp( const gsl_interp* interpolator, 
                          const Real* xTable, const Real* yTable,
                          const Real x )
{
    const Real y( gsl_interp_eval( interpolator, xTable, yTable, x, 0 ) );
    
    return y;
}


static void 
fillInterpolatorVector( std::vector<gsl_interp*>& interpolatorVector,
                        const sb_table::Table* const table[],
                        const UnsignedInteger N,
                        const gsl_interp_type* interpType )
{
    interpolatorVector.resize( N + 1 );
    for( UnsignedInteger n( 0 ); n <= N; ++n )
    {
        const sb_table::Table* tablen( table[n] ); 
        gsl_interp* interp = gsl_interp_alloc( interpType,
                                               tablen->N );
        
        gsl_interp_init( interp, tablen->x, tablen->y,
                         tablen->N );
        
        interpolatorVector[n] = interp;
        
        //printf("n i %d %d\n",n,tablen.size());
    }
}


const UnsignedInteger
SphericalBesselGenerator::getMaxNJ()
{
    return sb_table::sj_table_max;
}

const UnsignedInteger
SphericalBesselGenerator::getMaxNY()
{
    return sb_table::sy_table_max;
}

static const sb_table::Table* 
getSJTable( const UnsignedInteger n )
{
    return sb_table::sj_table[n];
}


static const sb_table::Table* 
getSYTable( const UnsignedInteger n )
{
    return sb_table::sy_table[n];
}

inline const Real SphericalBesselGenerator::_j_table( const UnsignedInteger n,
                                                      const Real z ) const
{
    const sb_table::Table* tablen( getSJTable( n ) );

    return interp( this->sjInterpolatorVector[n],
                   tablen->x, tablen->y, z );
}

inline const Real SphericalBesselGenerator::_y_table( const UnsignedInteger n, 
                                                      const Real z ) const
{
    const sb_table::Table* tablen( getSYTable( n ) );

    return interp( this->syInterpolatorVector[n],
                   tablen->x, tablen->y, z );
}


const Real 
SphericalBesselGenerator::j( const UnsignedInteger n, const Real z ) const
{
    if( n > getMaxNJ() )
    {
        return this->_j( n, z );
    }
    
    const sb_table::Table* table( getSJTable( n ) );
    
    const Real minz( table->x[2] );
    const Real maxz( table->x[table->N - 3] );
    
    if( z < minz )
    {
        return this->_j( n, z );
    }
    else if( z < maxz )
    {
        // table
        return this->_j_table( n, z );
    }
    else
    {
        return this->_j( n, z );
    }
}

const Real 
SphericalBesselGenerator::y( const UnsignedInteger n, const Real z ) const
{
    if( n > getMaxNY() )
    {
        return this->_y( n, z );
    }
    
    const sb_table::Table* table( getSYTable( n ) );
    
    const Real minz( table->x[0] );
    const Real maxz( table->x[table->N - 1] );
    
    if( z < minz )
    {
        return this->_y( n, z );
    }
    else if( z < maxz )
    {
        // table
        return this->_y_table( n, z );
    }
    else
    {
        return this->_y( n, z );
    }
}

void SphericalBesselGenerator::fillTables()
{
    fillInterpolatorVector( this->sjInterpolatorVector, sb_table::sj_table,
                            sb_table::sj_table_max, this->interpType );
    fillInterpolatorVector( this->syInterpolatorVector, sb_table::sy_table,
                            sb_table::sy_table_max, this->interpType );
}

