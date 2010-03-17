#ifndef REGION_HPP
#define REGION_HPP

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include <sstream>
#include "Structure.hpp"

template<typename Tid_, typename Tshape_>
class Region: public Structure<Tid_>
{
public:
    typedef Structure<Tid_> base_type;
    typedef Tid_ identifier_type;
    typedef Tshape_ shape_type;

public:
    virtual ~Region() {}

    identifier_type const& id() const
    {
        return id_;
    }

    identifier_type& id()
    {
        return id_;
    }

    shape_type& shape()
    {
        return shape_;
    }

    shape_type const& shape() const
    {
        return shape_;
    }

    bool operator==(Region const& rhs) const
    {
        return id_ == rhs.id() && shape_ == rhs.shape();
    }

    bool operator!=(Region const& rhs) const
    {
        return !operator==(rhs);
    }

    virtual std::size_t hash() const
    {
#if defined(HAVE_TR1_FUNCTIONAL)
        using std::tr1::hash;
#elif defined(HAVE_STD_HASH)
        using std::hash;
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
        using boost::hash;
#endif
        return hash<identifier_type>()(id()) ^ hash<shape_type>()(shape());
    }

    virtual std::string as_string() const
    {
        std::ostringstream out;
        out << "Region(" << id() << ":" << shape() << ")";
        return out.str(); }

    Region(identifier_type const& id, shape_type const& shape)
        : base_type(id), shape_(shape) {}

protected:
    identifier_type id_;
    shape_type shape_;
};

#endif /* REGION_HPP */
