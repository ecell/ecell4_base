#ifndef ECELL4_CORE_OBJECT_ID_CONTAINER_HPP
#define ECELL4_CORE_OBJECT_ID_CONTAINER_HPP
#include <ecell4/core/type_name_of.hpp>
#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <functional>
#include <utility>
#include <cassert>

namespace ecell4
{

//
// A container that combines map<id, index> and vector<{id, object}>.
// It does not check collision or any other kind of object dependent tests,
// but encapsulates relationships between id-idx map and a vector.
// It also has a SerialIDGenerator, but it does not automatically generate
// a new ID. World or some other higher order container does it.
//
template<typename Tid, typename Tobject>
class ObjectIDContainer
{
public:
    using identifier_type   = Tid;
    using object_type       = Tobject;
    using self_type         = ObjectIDContainer<identifier_type, object_type>;
    using id_idx_map_type   = std::unordered_map<identifier_type, std::size_t>;
    using value_type        = std::pair<identifier_type, object_type>;
    using container_type    = std::vector<value_type>;
    using size_type         = typename container_type::size_type;
    using difference_type   = typename container_type::difference_type;
    using iterator          = typename container_type::iterator;
    using const_iterator    = typename container_type::const_iterator;
    using id_generator_type = SerialIDGenerator<identifier_type>;

public:

    ObjectIDContainer()  = default;
    ~ObjectIDContainer() = default;
    ObjectIDContainer(const ObjectIDContainer&) = default;
    ObjectIDContainer(ObjectIDContainer&&)      = default;
    ObjectIDContainer& operator=(const ObjectIDContainer&) = default;
    ObjectIDContainer& operator=(ObjectIDContainer&&)      = default;

    bool update(const identifier_type& id, const object_type& obj)
    {
        const auto found(this->idxmap_.find(id));
        if(found == this->idxmap_.end())
        {
            idxmap_[id] = this->objects_.size();
            objects_.emplace_back(id, obj);
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
        const auto found(this->idxmap_.find(id));
        if(found == this->idxmap_.end())
        {
            throw_exceptions<NotFound>(utils::type_name_of<self_type>::value(),
                    "::remove(id=", id, "): object not found");
        }

        const size_type idx(found->second), last(objects_.size()-1);
        if(idx != last)
        {
            idxmap_[objects_[last].first] = idx;
            objects_[idx] = std::move(objects_[last]);
        }
        objects_.pop_back();
        idxmap_.erase(id);
        return;
    }

    identifier_type gen_id() {return this->idgen_();}

    std::pair<identifier_type, object_type> get(const identifier_type& id) const
    {
        const auto found(this->idxmap_.find(id));
        if(found == this->idxmap_.end())
        {
            throw_exceptions<NotFound>(utils::type_name_of<self_type>::value(),
                    "::get(id=", id, "): object not found");
        }
        return objects_[found->second];
    }

    bool has(const identifier_type& id) const {return idxmap_.count(id) != 0;}
    container_type const& list() const noexcept {return objects_;}

    std::size_t size() const noexcept {return objects_.size();}
    bool       empty() const noexcept {return objects_.empty();}

    bool diagnosis() const
    {
        for(const auto& item : this->idxmap_);
        {
            const auto& id  = item.first;
            const auto& idx = item.second;

            identifier_type id_by_cont;
            try
            {
                id_by_cont = objects_.at(idx).first;
            }
            catch(std::out_of_range& oor)
            {
                throw_exceptions<IllegalState>(utils::type_name_of<self_type>::value(),
                    "::diagnosis: object(id=", id, ") not found in the container");
            }
            if(id_by_cont != id)
            {
                throw_exceptions<IllegalState>(utils::type_name_of<self_type>::value(),
                    "::diagnosis: object(id=", id, ") index is invalid: ", idx,
                    "-th object has id ", id_by_cont);
            }
        }
        for(std::size_t i=0; i<objects_.size(); ++i)
        {
            const auto& id = objects_.at(i).first;
            if(this->idxmap_.count(id) == 0)
            {
                throw_exceptions<IllegalState>(utils::type_name_of<self_type>::value(),
                    "::diagnosis: object(id=", id, ") not found in idxmap");
            }
            if(this->idxmap_.at(id) != i)
            {
                throw_exceptions<IllegalState>(utils::type_name_of<self_type>::value(),
                    "::diagnosis: object(id=", id, ") has different idx (", idxmap_.at(id),
                    ") in idxmap");
            }
        }
        return true;
    }

    value_type& at(const identifier_type& id)
    {
        const auto found(idxmap_.find(id));
        if(found == idxmap_.end())
        {
            throw_exceptions<NotFound>(utils::type_name_of<self_type>::value(),
                    "::at(id=", id, "): object not found");
        }
        return objects_.at(found->second);
    }
    value_type const& at(const identifier_type& id) const
    {
        const auto found(idxmap_.find(id));
        if(found == idxmap_.end())
        {
            throw_exceptions<NotFound>(utils::type_name_of<self_type>::value(),
                    "::at(id=", id, "): object not found");
        }
        return objects_.at(found->second);
    }

    value_type& operator[](const identifier_type& id)
    {
        return objects_[idxmap_[id]];
    }
    value_type const& operator[](const identifier_type& id) const
    {
        return objects_[idxmap_[id]];
    }

    iterator       begin()        noexcept {return objects_.begin();}
    iterator       end()          noexcept {return objects_.end();}
    const_iterator begin()  const noexcept {return objects_.begin();}
    const_iterator end()    const noexcept {return objects_.end();}
    const_iterator cbegin() const noexcept {return objects_.begin();}
    const_iterator cend()   const noexcept {return objects_.end();}

private:

    id_generator_type idgen_;
    id_index_map_type idxmap_;
    container_type    objects_;
};
} // ecell4
#endif// ECELL4_OBJECT_ID_CONTAINER
