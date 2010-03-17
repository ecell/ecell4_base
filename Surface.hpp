#ifndef SURFACE_HPP
#define SURFACE_HPP

#include "Vector3.hpp"
#include "Single.hpp"
#include "Pair.hpp"
#include "Multi.hpp"

template<typename Tid_, typename Tpos_, typename Tlen_ = typename element_type_of<Tpos_>::type>
class Surface
{
public:
    typedef Tid_ identifier_type;
    typedef Tpos_ position_type;
    typedef Tlen_ length_type;

public:
    virtual ~Surface() {}

    identifier_type const& id() const
    {
        return id_;
    }

    identifier_type& id()
    {
        return id_;
    }

    Surface(identifier_type const& id): id_(id) {}

protected:
    identifier_type id_;
};

#endif /* SURFACE_HPP */
