#ifndef __ECELL4_SHAPE_CONTAINER_HPP
#define __ECELL4_SHAPE_CONTAINER_HPP

#include "config.h"

#include <ostream>
#if defined(HAVE_TR1_FUNCTIONAL)
#include <tr1/functional>
#elif defined(HAVE_STD_HASH)
#include <functional>
#elif defined(HAVE_BOOST_FUNCTIONAL_HASH_HPP)
#include <boost/functional/hash.hpp>
#endif


#include <utility>

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "Real3.hpp"
#include "Identifier.hpp"
#include "Shape.hpp"

#include <boost/shared_ptr.hpp>

#include <vector>
#include <map>

namespace ecell4
{

struct ShapeID:
    public Identifier<ShapeID, unsigned long long, int>
{
    typedef Identifier<ShapeID, unsigned long long, int> base_type;
    ShapeID(const value_type& value = value_type(0, 0))
        :base_type(value)
    {
        ;
    }
};

class ShapeContainer 
{
public:
    typedef std::vector<std::pair<ShapeID, boost::shared_ptr<Shape> > >
        shape_container_type;

protected:
    typedef utils::get_mapper_mf<
        ShapeID, shape_container_type::size_type>::type shape_map_type;
public:
    // Constructor
    ShapeContainer(void);

    Integer num_shape(void) const
    {
        return shapes_.size();
    }
    Integer num_shape(Shape::dimension_kind &dim) const
    {
        Integer ret = 0;
        for(shape_container_type::const_iterator it = shapes_.begin(); it != shapes_.end(); it++)
        {
            if (it->second->dimension() == dim) {
                ret++;
            }
        }
        return ret;
    }

    std::pair<ShapeID, boost::shared_ptr<Shape> > get_shape(const ShapeID& shape_id) const
    {
        shape_map_type::const_iterator i(index_map_.find(shape_id));
        if (i != index_map_.end())
        {
            throw NotFound("shape not found");
        }
        return shapes_[(*i).second];
    }

    shape_container_type list_shapes() const
    {
        return shapes_;
    }

    shape_container_type list_shapes(Shape::dimension_kind const &dim) const 
    {
        shape_container_type ret;
        for(shape_container_type::const_iterator it = shapes_.begin(); it != shapes_.end(); it++)
        {
            if (it->second->dimension() == dim) {
                ret.push_back(*it);
            }
        }
        return ret;
    }
    bool has_shape(ShapeID const &shape_id) const
    {
        shape_map_type::const_iterator i(index_map_.find(shape_id));
        return (i != index_map_.end());
    }

    bool update_shape(const ShapeID& shape_id, const boost::shared_ptr<Shape> s)
    {
        shape_map_type::const_iterator i(index_map_.find(shape_id));
        if (i == index_map_.end())
        {
            shape_container_type::size_type idx(shapes_.size());
            index_map_[shape_id] = idx;
            shapes_.push_back(std::make_pair(shape_id, s));
            return true;
        }
        else
        {
            shapes_[(*i).second] = std::make_pair(shape_id, s);
            return false;
        }
    }
protected:
    shape_container_type shapes_;
    shape_map_type index_map_;
};

}   //ecell4

#endif
