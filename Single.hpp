#ifndef SINGLE_HPP
#define SINGLE_HPP

#include "Domain.hpp"

template<typename T_>
class Single: public Domain<T_>
{
public:
    typedef Domain<T_> base_type;
};

#endif /* SINGLE_HPP */
