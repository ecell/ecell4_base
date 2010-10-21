#ifndef SURFACE_HPP
#define SURFACE_HPP

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif

#include <sstream>
#include "ParticleSimulationStructure.hpp"
#include "Box.hpp"
#include "Cylinder.hpp"
#include "Sphere.hpp"
#include "Plane.hpp"

template<typename Ttraits_>
class Surface: public ParticleSimulationStructure<Ttraits_>
{
public:
    typedef ParticleSimulationStructure<Ttraits_> base_type;
    typedef typename base_type::identifier_type identifier_type;
    typedef typename base_type::length_type length_type;

public:
    virtual ~Surface() {}

    Surface(identifier_type const& id): base_type(id) {}

    virtual length_type minimal_distance(length_type const& radius) const = 0;
};

template<typename Ttraits_, typename Tshape_>
class BasicSurfaceImpl: public Surface<Ttraits_>
{
public:
    typedef Surface<Ttraits_> base_type;
    typedef Tshape_ shape_type;
    typedef typename base_type::identifier_type identifier_type;
    typedef typename base_type::length_type length_type;
    typedef typename base_type::position_type position_type;

public:
    virtual ~BasicSurfaceImpl() {}

    shape_type& shape()
    {
        return shape_;
    }

    shape_type const& shape() const
    {
        return shape_;
    }

    virtual bool operator==(Structure<typename Ttraits_::world_type::traits_type> const& rhs) const
    {
        BasicSurfaceImpl const* _rhs(dynamic_cast<BasicSurfaceImpl const*>(&rhs));
        return _rhs && base_type::id_ == rhs.id() && shape_ == _rhs->shape();
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
        return hash<identifier_type>()(base_type::id_) ^ hash<shape_type>()(shape());
    }

    virtual std::string as_string() const
    {
        std::ostringstream out;
        out << "Surface(" << base_type::id_ << ":" << shape() << ")";
        return out.str();
    }

    //std::pair<position_type, length_type>
    position_type projected_point(position_type const& pos) const
    {
        // Todo. Return pair and Pythonify.
        return ::projected_point(shape(), pos).first;
    }

    BasicSurfaceImpl(identifier_type const& id, shape_type const& shape)
        : base_type(id), shape_(shape) {}

protected:
    shape_type shape_;
};

#endif /* SURFACE_HPP */
