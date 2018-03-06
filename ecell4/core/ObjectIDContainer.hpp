#ifndef ECELL4_OBJECT_ID_CONTAINER
#define ECELL4_OBJECT_ID_CONTAINER
#include <ecell4/core/get_mapper_mf.hpp>
#include <ecell4/core/type_name_of.hpp>
#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <boost/format.hpp>
#include <functional>
#include <utility>
#include <cassert>

namespace ecell4
{

//! @brief utility container to contain an object (e.g. Particle) with its ID.
template<typename T_id, typename T_obj>
class ObjectIDContainer
{
public:
    typedef T_id    identifier_type;
    typedef T_obj   object_type;
    typedef ObjectIDContainer<T_id, T_obj> self_type;
    typedef typename utils::get_mapper_mf<identifier_type, std::size_t>::type
            id_index_map_type;

    typedef std::pair<identifier_type, object_type>  value_type;
    typedef std::vector<value_type>                  container_type;
    typedef typename container_type::size_type       size_type;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::iterator        iterator;
    typedef typename container_type::const_iterator  const_iterator;

public:

    bool update(const identifier_type& id, const object_type& obj)
    {
        const typename id_index_map_type::iterator found(this->idxmap_.find(id));
        if(found == this->idxmap_.end())
        {
            const size_type idx(this->objects_.size());
            idxmap_[id] = idx;
            objects_.push_back(std::make_pair(id, obj));
            return true;
        }
        else
        {
            objects_[found->second] = std::make_pair(id, obj);
            return false;
        }
    }
    void remove(const identifier_type& id)
    {
        const typename id_index_map_type::iterator found(this->idxmap_.find(id));
        if(found == this->idxmap_.end())
        {
            throw NotFound((boost::format(
                    "%1%::remove: object(id=%2%) not found") %
                    (utils::type_name_of<self_type>::value()) %
                    id).str());
        }
        const size_type idx(found->second), last(objects_.size()-1);
        if(idx != last)
        {
            idxmap_[objects_[last].first] = idx;
            objects_[idx] = objects_[last];
        }
        objects_.pop_back();
        idxmap_.erase(id);
    }

    bool has(const identifier_type& id) const {return idxmap_.count(id) == 1;}
    container_type const& list() const throw() {return objects_;}

    std::size_t size() const throw() {return objects_.size();}
    bool       empty() const throw() {return objects_.empty();}

    bool diagnosis() const
    {
        for(typename id_index_map_type::const_iterator
                i(idxmap_.begin()), e(idxmap_.end()); i!=e; ++i)
        {
            identifier_type id_by_cont;
            try
            {
                id_by_cont = objects_.at(i->second);
            }
            catch(std::out_of_range& oor)
            {
                throw std::runtime_error((boost::format(
                    "%1%::diagnosis: object(id=%2%) not found in the container") %
                    (utils::type_name_of<self_type>::value()) %
                    i->first).str());
            }
            if(id_by_cont != i->first)
            {
                 throw std::runtime_error((boost::format(
                    "%1%::diagnosis: object(id=%2%) index invalid: "
                    "index %3% has different id %4%.") %
                    (utils::type_name_of<self_type>::value()) %
                    i->first % i->second % id_by_cont).str());
            }
        }

        for(std::size_t i=0, e(this->objects_.size()); i<e; ++i)
        {
            if(this->idxmap_.count(this->objects_.at(i).first) == 0)
            {
                throw std::runtime_error((boost::format(
                    "%1%::diagnosis: object(id=%2%) not found in the map") %
                    (utils::type_name_of<self_type>::value()) %
                    this->objects_.at(i).first).str());
            }
        }
        return true;
    }

    value_type& at(const identifier_type& id)
    {
        typename id_index_map_type::const_iterator found(idxmap_.find(id));
        if(found == idxmap_.end())
        {
            throw std::out_of_range((
                boost::format("%1%::at: no element with id (%2%)") %
                (utils::type_name_of<self_type>::value()) % id).str());
        }
        return objects_[found->second];
    }
    value_type const& at(const identifier_type& id) const
    {
        // C++98 std::map does not have std::map::at() const.
        typename id_index_map_type::const_iterator found(idxmap_.find(id));
        if(found == idxmap_.end())
        {
            throw std::out_of_range((
                boost::format("%1%::at const: no element with id (%2%)") %
                (utils::type_name_of<self_type>::value()) % id).str());
        }
        return objects_.at(found->second);
    }

    value_type& operator[](const identifier_type& id)
    {
        return objects_[idxmap_.at(id)];
    }
    value_type const& operator[](const identifier_type& id) const
    {
        return objects_[idxmap_.at(id)];
    }

    value_type&       front()       throw() {return objects_.front();}
    value_type const& front() const throw() {return objects_.front();}
    value_type&       back()        throw() {return objects_.back();}
    value_type const& back()  const throw() {return objects_.back();}

    iterator       begin()        throw() {return objects_.begin();}
    iterator       end()          throw() {return objects_.end();}
    const_iterator begin()  const throw() {return objects_.begin();}
    const_iterator end()    const throw() {return objects_.end();}
    const_iterator cbegin() const throw() {return objects_.begin();}
    const_iterator cend()   const throw() {return objects_.end();}

private:

    id_index_map_type idxmap_;
    container_type    objects_;
};
} // ecell4
#endif// ECELL4_OBJECT_ID_CONTAINER
