#ifndef UTILS_HPP
#define UTILS_HPP

#include <map>
#include <vector>

namespace get_default_impl
{
    namespace std
    {
        template<typename Tkey_, typename Tval_>
        struct map
        {
            typedef ::std::map<Tkey_, Tval_> type;
        };

        template<typename Tval_>
        struct vector
        {
            typedef ::std::vector<Tval_> type; 
        };
    } // std
} // namespace get_default_impl

void gsl_error_handler( char const* reason, char const* file, int line, int gsl_errno );

#endif /* UTILS_HPP */
