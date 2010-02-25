#ifndef UNASSIGNABLE_ADAPTER_HPP
#define UNASSIGNABLE_ADAPTER_HPP

#include <algorithm>
#include <boost/utility/enable_if.hpp>
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
#include <boost/iterator/transform_iterator.hpp>
#include "utils/fun_wrappers.hpp"

template<typename T_, template<typename> class TT_>
struct unassignable_adapter
{
private:
    struct placeholder { char _[sizeof(T_)]; };

public:
    typedef typename TT_<placeholder>::type container_type;
    typedef T_ value_type;

    typedef value_type* pointer;
    typedef value_type const* const_pointer;
    typedef value_type& reference;
    typedef value_type const& const_reference;

private:
    typedef reinterpret_caster<reference, typename container_type::value_type&> caster;
    typedef reinterpret_caster<const_reference, typename container_type::value_type const&> const_caster;

public:

    typedef typename boost::range_size<container_type>::type size_type;
    typedef typename boost::range_difference<container_type>::type difference_type;
    typedef typename boost::transform_iterator<caster,
        typename boost::range_iterator<container_type>::type> iterator;
    typedef typename boost::transform_iterator<const_caster,
        typename boost::range_const_iterator<container_type>::type> const_iterator;
    typedef typename boost::transform_iterator<caster,
        typename boost::range_reverse_iterator<container_type>::type> reverse_iterator;
    typedef typename boost::transform_iterator<const_caster,
        typename boost::range_const_reverse_iterator<container_type>::type> const_reverse_iterator;

    size_type size() const
    {
        return boost::size(cntnr_);
    }

    iterator begin()
    {
        return iterator(boost::begin(cntnr_), caster());
    }

    const_iterator begin() const
    {
        return const_iterator(boost::begin(cntnr_), const_caster());
    }

    iterator end()
    {
        return iterator(boost::end(cntnr_), caster());
    }

    const_iterator end() const
    {
        return const_iterator(boost::end(cntnr_), const_caster());
    }

    reverse_iterator rbegin()
    {
        return reverse_iterator(boost::rbegin(cntnr_), caster());
    }

    const_reverse_iterator rbegin() const
    {
        return const_reverse_iterator(boost::rbegin(cntnr_), caster());
    }

    reverse_iterator rend()
    {
        return reverse_iterator(boost::rend(cntnr_), caster());
    }

    const_reverse_iterator rend() const
    {
        return const_reverse_iterator(boost::rend(cntnr_), const_caster());
    }

    void clear()
    {
        cntnr_.clear();
    }

    void resize(size_type n, T_ const& v)
    {
        cntnr_.resize(n, reinterpret_cast<typename container_type::value_type const&>(v));
    }

    value_type const& at(size_type const& pos) const
    {
        return reinterpret_cast<value_type const&>(cntnr_.at(pos));
    }

    value_type& at(size_type const& pos)
    {
        return reinterpret_cast<value_type&>(cntnr_.at(pos));
    }

    value_type const& operator[](size_type const& pos) const
    {
        return reinterpret_cast<value_type const&>(cntnr_[pos]);
    }

    void set(size_type const& pos, value_type const& v)
    {
        cntnr_[pos] = reinterpret_cast<typename container_type::value_type const&>(v);
    }

    void insert(iterator const& pos, value_type const& v)
    {
        cntnr_.insert(pos, reinterpret_cast<typename container_type::value_type const&>(v));
    }

    void insert(iterator const& pos, size_type n, value_type const v)
    {
        cntnr_.insert(pos, n, reinterpret_cast<typename container_type::value_type const&>(v));
    }

    template<typename Titer_>
    void insert(iterator const& pos, Titer_ const& b, Titer_ const& e)
    {
        typedef reinterpret_caster<typename container_type::value_type const&, const_reference> reverse_caster;
        typedef boost::transform_iterator<reverse_caster, Titer_> transform_iterator;
        cntnr_.insert(pos, transform_iterator(b, reverse_caster()), transform_iterator(e, caster()));
    }

    void push_front(value_type const& v)
    {
        cntnr_.push_front(reinterpret_cast<typename container_type::value_type const&>(v));
    }

    void push_back(value_type const& v)
    {
        cntnr_.push_back(reinterpret_cast<typename container_type::value_type const&>(v));
    }

    void pop_front()
    {
        cntnr_.pop_front();
    }

    void pop_back()
    {
        cntnr_.pop_back();
    }

    reference front()
    {
        return reinterpret_cast<reference>(cntnr_.front());
    }

    const_reference front() const
    {
        return reinterpret_cast<const_reference>(cntnr_.front());
    }

    reference back()
    {
        return reinterpret_cast<reference>(cntnr_.back());
    }

    const_reference back() const
    {
        return reinterpret_cast<const_reference>(cntnr_.back());
    }

    iterator erase(iterator const& pos)
    {
        return cntnr_.erase(pos);
    }

    iterator erase(iterator const& b, iterator const& e)
    {
        return cntnr_.erase(b, e);
    }

    unassignable_adapter() {}

private:
    container_type cntnr_;
};

#endif /* UNASSIGNABLE_ADAPTER_HPP */
