#ifndef ECELL4_DYNAMICPRIORITYQUEUE_HPP
#define ECELL4_DYNAMICPRIORITYQUEUE_HPP
//
// written by Koichi Takahashi based on the initial version by Eiichiro Adachi.
// modified by Mozoyoshi Koizumi
//

#include <ecell4/core/config.h>

#include <functional>
#include <vector>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <cstring>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

//#define HAVE_TR1_UNORDERED_MAP

#if HAVE_UNORDERED_MAP
#include <unordered_map>
#elif HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
#elif HAVE_BOOST_UNORDERED_MAP_HPP
#include <boost/unordered_map.hpp>
#else
#include <map>
#endif /* HAVE_UNORDERED_MAP */

#ifdef DEBUG
#include <iostream>
#endif

#include "swap.hpp"


namespace ecell4
{

template<typename Tid_>
struct default_id_generator
{
    typedef Tid_ identifier_type;

    default_id_generator(): next_() {}

    default_id_generator(identifier_type const& first): next_(first) {}

    identifier_type operator()()
    {
        return ++next_;
    }

protected:
    identifier_type next_;
};


template<typename Tid_ = unsigned long long,
         typename Tindex_ = std::size_t,
         typename Tidgen_ = default_id_generator<Tid_> >
class persistent_id_policy
{
public:
    typedef Tid_ identifier_type;
    typedef Tindex_ index_type;
    typedef Tidgen_ identifier_generator;

protected:
    struct hasher
        : public std::unary_function<identifier_type, std::size_t>
    {
        std::size_t operator()(identifier_type value) const
        {
            return static_cast<std::size_t>(value) ^
                static_cast<std::size_t>(
                    value >> (sizeof(identifier_type) * 8 / 2));
        }
    };
#if HAVE_UNORDERED_MAP
    typedef std::unordered_map<identifier_type, index_type, hasher> index_map;
#elif HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<identifier_type, index_type, hasher> index_map;
#elif HAVE_BOOST_UNORDERED_MAP_HPP
    typedef boost::unordered_map<identifier_type, index_type, hasher> index_map;
#else
    typedef std::map<identifier_type, index_type> index_map;
#endif

public:
    index_type index(identifier_type const& id) const
    {
        typename index_map::const_iterator i(index_map_.find(id));
        if (i == index_map_.end())
        {
            throw std::out_of_range((boost::format("%s: Key not found (%s)")
                % __FUNCTION__ % boost::lexical_cast<std::string>(id)).str());
            // throw std::out_of_range((boost::format("%s: Key not found (%s)")
            //     % __PRETTY_FUNCTION__ % boost::lexical_cast<std::string>(id)).str());
        }
        return (*i).second;
    }

    identifier_type push(index_type index)
    {
        const identifier_type id(idgen_());
        index_map_.insert(typename index_map::value_type(id, index));
        return id;
    }

    void pop(index_type index, identifier_type id, identifier_type last_item_id)
    {
        index_map_[last_item_id] = index;
        index_map_.erase(id);
    }

    void clear()
    {
        index_map_.clear();
    }

private:
    index_map index_map_;
    identifier_generator idgen_;
};

template<typename Tindex_ = std::size_t>
class volatile_id_policy
{
public:
    typedef Tindex_ identifier_type;
    typedef Tindex_ index_type;

    index_type index(identifier_type const& id) const
    {
        return id;
    }

    identifier_type push(index_type index)
    {
        return index;
    }

    void pop(index_type, identifier_type, identifier_type) {}

    void clear() {}
};


/**
   Dynamic priority queue for items of type Titem_.

   When Tpolicy_ template parameter is persistent_id_policy, identifier_types assigned
   to pushed items are persistent for the life time of this priority
   queue.

   When Volatileidentifier_typePolicy template parameter is used as the Tpolicy_,
   identifier_types are valid only until the next call of pop or push methods.
   However, Volatileidentifier_typePolicy saves some memory and eliminates the
   overhead incurred in pop/push methods.
*/

template<typename Titem_, typename Tcomparator = std::less_equal<Titem_>, class Tpolicy_ = persistent_id_policy<> >
class DynamicPriorityQueue: private Tpolicy_
{
public:
    typedef Tpolicy_ policy_type;
    typedef typename policy_type::identifier_type identifier_type;
    typedef typename policy_type::index_type index_type;
    typedef Titem_ element_type;
    typedef std::pair<identifier_type, element_type> value_type;
    typedef Tcomparator comparator_type;

protected:
    typedef std::vector<value_type> value_vector;
    typedef std::vector<index_type> index_vector;

public:
    typedef typename value_vector::size_type size_type;
    typedef typename value_vector::const_iterator iterator;
    typedef typename value_vector::const_iterator const_iterator;

public:
    bool empty() const
    {
        return items_.empty();
    }

    size_type size() const
    {
        return items_.size();
    }

    void clear();

    value_type const& top() const
    {
        return items_[top_index()];
    }

    value_type const& second() const
    {
        return items_[second_index()];
    }

    element_type const& get(identifier_type id) const
    {
        return items_[policy_type::index(id)].second;
    }

    void pop()
    {
        pop_by_index(top_index());
    }

    void pop(identifier_type id)
    {
        pop_by_index(policy_type::index(id));
    }

    void replace(value_type const& item);

    identifier_type push(element_type const& item);

    element_type const& operator[](identifier_type id) const
    {
        return get(id);
    }

    const_iterator begin() const
    {
        return items_.begin();
    }

    const_iterator end() const
    {
        return items_.end();
    }

    // self-diagnostic methods
    bool check() const; // check all
    bool check_size() const;
    bool check_position_mapping() const;
    bool check_heap() const;


protected:
    index_type top_index() const
    {
        return heap_[0];
    }

    index_type second_index() const
    {
        if (size() <= 1)
        {
            throw std::out_of_range("DynamicPriorityQueue::second_index():"
                                     " item count less than 2.");
        }

        const index_type index1(heap_[1]);

        if (size() == 2)
        {
            return index1;
        }

        const index_type index2(heap_[2]);
        if (comp(items_[index1].second, items_[index2].second))
        {
            return index1;
        }
        else
        {
            return index2;
        }
    }

    void pop_by_index(index_type index);

    void move(index_type index)
    {
        const index_type pos(position_vector_[index]);
        move_pos(pos);
    }

    void move_top()
    {
        move_down_pos(0);
    }

    void move_pos(index_type pos);

    void move_up(index_type index)
    {
        const index_type position(position_vector_[index]);
        move_up_pos(position);
    }

    void move_down(index_type index)
    {
        const index_type position(position_vector_[index]);
        move_down_pos(position);
    }


    void move_up_pos(index_type position, index_type start = 0);
    void move_down_pos(index_type position);

    void move_up_pos_impl(index_type position, index_type start = 0);
    void move_down_pos_impl(index_type position);

private:
    value_vector items_;
    index_vector heap_;
    index_vector position_vector_;

    comparator_type comp;
};

template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline void DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::clear()
{
    items_.clear();
    heap_.clear();
    position_vector_.clear();
    policy_type::clear();
}


template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline void DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::move_pos(index_type pos)
{
    const index_type index(heap_[pos]);
    const value_type& item(items_[index]);
    const index_type succ(2 * pos + 1);
    if (succ < size())
    {
        if (comp(items_[heap_[succ]].second, item.second) || (succ + 1 < size() && comp(items_[heap_[succ + 1]].second, item.second)))
        {
            move_down_pos_impl(pos);
            return;
        }
    }

    move_up_pos(pos);
}

template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline void DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::move_up_pos(index_type position, index_type start)
{
    if (position == 0)
        return;

    const index_type index(heap_[position]);
    const value_type& item(items_[index]);

    const index_type pred((position - 1) / 2);
    const index_type predindex_type(heap_[pred]);

    if (comp(item.second, items_[predindex_type].second))
    {
        move_up_pos_impl(position, start);
    }
}


template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline void DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::move_down_pos(index_type position)
{
    const index_type index(heap_[position]);
    const value_type& item(items_[index]);

    const index_type succ(2 * position + 1);
    if (succ < size())
    {
        if (comp(items_[heap_[succ]].second, item.second) || (succ + 1 < size() && comp(items_[heap_[succ + 1]].second, item.second)))
        {
            move_down_pos_impl(position);
        }
    }
}

template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline void DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::move_up_pos_impl(index_type position, index_type start)
{
    const index_type index(heap_[position]);
    const value_type& item(items_[index]);

    if (position <= start)
    {
        return;
    }

    index_type pos(position);
    index_type pred((pos - 1) / 2);
    index_type predindex_type(heap_[pred]);

    do
    {
        heap_[pos] = predindex_type;
        position_vector_[predindex_type] = pos;
        pos = pred;

        if (pos <= start)
        {
            break;
        }

        pred = (pos - 1) / 2;
        predindex_type = heap_[pred];

    } while (! comp(items_[predindex_type].second, item.second));

    heap_[pos] = index;
    position_vector_[index] = pos;
}


template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline void DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::move_down_pos_impl(index_type position)
{
    const index_type index(heap_[position]);

    index_type succ(2 * position + 1);
    index_type pos(position);
    while (succ < size())
    {
        const index_type right_pos(succ + 1);
        if (right_pos < size() && !comp(items_[heap_[succ]].second, items_[heap_[right_pos]].second))
        {
            succ = right_pos;
        }

        heap_[pos] = heap_[succ];
        position_vector_[heap_[pos]] = pos;
        pos = succ;
        succ = 2 * pos + 1;
    }

    heap_[pos] = index;
    position_vector_[index] = pos;

    move_up_pos(pos, position);
}



template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline typename DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::identifier_type
DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::push(Titem_ const& item)
{
    const index_type index(items_.size());
    const identifier_type id(policy_type::push(index));
    items_.push_back(value_type(id, item));
    // index == pos at this time.
    heap_.push_back(index);
    position_vector_.push_back(index);
    move_up_pos(index);
    return id;
}


template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline void DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::pop_by_index(index_type index)
{
    value_type& item(items_[index]);
    // 1. update index<->identifier_type mapping.
    policy_type::pop(index, item.first, items_.back().first);

    // 2. pop the item from the items_.
    blit_swap(item, items_.back());
    items_.pop_back();

    const index_type removed_pos(position_vector_[index]);
    const index_type moved_pos(position_vector_.back());

    // 3. swap position_vector_[end] and position_vector_[index]
    position_vector_[index] = moved_pos;
    heap_[moved_pos] = index;

    // 4. if heap_[end] and heap_[removed] do not overlap,
    //    swap these, pop back, and update the heap_.
    if (removed_pos != heap_.size() - 1)
    {
        heap_[removed_pos] = heap_.back();
        position_vector_[heap_.back()] = removed_pos;

        position_vector_.pop_back();
        heap_.pop_back();

        move_pos(removed_pos);
    }
    else  // if heap_[end] and heap_[removed] are the same, simply pop back.
    {
        position_vector_.pop_back();
        heap_.pop_back();
    }
}

template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline void DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::replace(value_type const& value)
{
    const index_type index(policy_type::index(value.first));
    items_[index].second = value.second;
    move(index);
}

template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline bool DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::check() const
{
    bool result(true);

    result = result && check_size();
    result = result && check_position_mapping();
    result = result && check_heap();

    return result;
}


template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline bool DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::check_size() const
{
    bool result(true);

    // check sizes of data structures.
    result = result && items_.size() == size();
    result = result && heap_.size() == size();
    result = result && position_vector_.size() == size();

    return result;
}


template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline bool DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::check_position_mapping() const
{
    bool result(true);

    // assert correct mapping between the heap_ and the position_vector_.
    for (index_type i(0); i < size(); ++i)
    {
        result = result && heap_[i] < size();
        result = result && position_vector_[i] < size();
        result = result && heap_[position_vector_[i]] == i;
    }

    return result;
}

template<typename Titem_, typename Tcomparator_, typename Tpolicy_>
inline bool DynamicPriorityQueue<Titem_, Tcomparator_, Tpolicy_>::check_heap() const
{
    bool result(true);

    // assert correct ordering of items in the heap_.

    for (index_type pos(0); pos < size(); ++pos)
    {
        const value_type& item(items_[heap_[pos]]);

        const index_type succ(pos * 2 + 1);
        if (succ < size())
        {
            result = result &&
                comp(item.second, items_[heap_[succ]].second);

            const index_type right_pos(succ + 1);
            if (right_pos < size())
            {
                result = result && comp(item.second, items_[heap_[right_pos]].second);
            }
        }

    }

    return result;
}

} // ecell4

#endif // ECELL4_DYNAMICPRIORITYQUEUE_HPP
