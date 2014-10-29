#ifndef MAP_ADAPTER_HPP
#define MAP_ADAPTER_HPP

#include <utility>
#include "utils.hpp"
#include <boost/range/size.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/rbegin.hpp>
#include <boost/range/rend.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/size_type.hpp>
#include <boost/range/difference_type.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/reverse_iterator.hpp>
#include <boost/range/const_reverse_iterator.hpp>
#include "assoc_container_traits.hpp"

template<typename Tcntnr_, typename Thdlr_>
struct map_adapter
{
public:
    typedef Tcntnr_ container_type;
    typedef typename boost::range_value<container_type>::type value_type;
    typedef typename boost::range_size<container_type>::type size_type;
    typedef typename boost::range_difference<container_type>::type difference_type;
    typedef typename boost::range_iterator<container_type>::type iterator;
    typedef typename boost::range_const_iterator<container_type>::type const_iterator;
    typedef typename boost::range_reverse_iterator<container_type>::type reverse_iterator;
    typedef typename boost::range_const_reverse_iterator<container_type>::type const_reverse_iterator;
    typedef typename assoc_key<Tcntnr_>::type key_type;
    typedef typename assoc_mapped<Tcntnr_>::type mapped_type;

    typedef value_type* pointer;
    typedef value_type const* const_pointer;
    typedef value_type& reference;
    typedef value_type const& const_reference;

public:
    ~map_adapter()
    {
        hdlr_.destroy(*this);
    }

    size_type size() const
    {
        return boost::size(cntnr_);
    }

    bool empty() const
    {
        return size() == 0;
    }

    iterator begin()
    {
        return boost::begin(cntnr_);
    }

    const_iterator begin() const
    {
        return boost::begin(cntnr_);
    }

    iterator end()
    {
        return boost::end(cntnr_);
    }

    const_iterator end() const
    {
        return boost::end(cntnr_);
    }

    reverse_iterator rbegin()
    {
        return boost::rbegin(cntnr_);
    }

    const_reverse_iterator rbegin() const
    {
        return boost::rbegin(cntnr_);
    }

    reverse_iterator rend()
    {
        return boost::rend(cntnr_);
    }

    const_reverse_iterator rend() const
    {
        return boost::end(cntnr_);
    }

    iterator find(key_type const& k)
    {
        return cntnr_.find(k);
    }

    const_iterator find(key_type const& k) const
    {
        return cntnr_.find(k);
    }

    iterator lower_bound(key_type const& k)
    {
        return cntnr_.lower_bound(k);
    }

    const_iterator lower_bound(key_type const& k) const
    {
        return cntnr_.lower_bound(k);
    }

    iterator upper_bound(key_type const& k)
    {
        return cntnr_.upper_bound(k);
    }

    const_iterator upper_bound(key_type const& k) const
    {
        return cntnr_.upper_bound(k);
    }

    size_type erase(key_type const& k)
    {
        return cntnr_.erase(k);
    }

    void erase(iterator const& pos)
    {
        return cntnr_.erase(pos);
    }

    void erase(iterator const& b, iterator const& e)
    {
        return cntnr_.erase(b, e);
    }

    std::pair<iterator, bool> insert(value_type const& v)
    {
        hdlr_.template insert<map_adapter>(v);
        return cntnr_.insert(v);
    }

    iterator insert(iterator const& hint, value_type const& v)
    {
        hdlr_.template insert<map_adapter>(v);
        return cntnr_.insert(hint, v);
    }

    template<typename Titer_>
    void insert(Titer_ const& b, Titer_ const& e)
    {
        hdlr_.template insert<map_adapter>(b, e);
        cntnr_.insert(b, e);
    }

    mapped_type const& at(key_type const& k) const
    {
        return cntnr_.at(k);
    }

    mapped_type& at(key_type const& k)
    {
        return cntnr_.at(k);
    }

    mapped_type const& operator[](key_type const& k) const
    {
        return at(k);
    }

    mapped_type& operator[](key_type& k) const
    {
        return at(k);
    }

    void clear()
    {
        cntnr_.clear();
    }

    map_adapter(Thdlr_ hdlr): hdlr_(hdlr) {}

private:
    container_type cntnr_;
    Thdlr_ hdlr_;
};

#endif /* MAP_ADAPTER_HPP */
