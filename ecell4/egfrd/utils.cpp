#include <stdexcept>
#include <string>
#include <boost/lexical_cast.hpp>
#include "utils.hpp"

// GSL error handler.
void gsl_error_handler( char const* reason, char const* file, int line, int gsl_errno )
{
    throw std::runtime_error( std::string( "GSL error: " ) +
                              std::string( reason ) +
                              std::string( " at " ) +
                              std::string( file ) + std::string( ":" ) +
                              boost::lexical_cast< std::string >( line ) );
}
