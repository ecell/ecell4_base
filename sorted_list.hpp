#ifndef _SORTED_LIST
#define _SORTED_LIST

#include <algorithm>
#include <functional>
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
#include "utils/fun_composition.hpp"

template<typename Tcntnr_, typename TweakOrdering_ = std::less<typename Tcntnr_::value_type> >
class sorted_list
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

    typedef value_type* pointer;
    typedef value_type const* const_pointer;
    typedef value_type& reference;
    typedef value_type const& const_reference;

    size_type size() const
    {
        return boost::size(cntnr_);
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

    void push(value_type const& v)
    {
        iterator i(std::upper_bound(cntnr_.begin(), cntnr_.end(), v,
                static_cast<TweakOrdering_ const&>(ord_)));
        cntnr_.insert(i, v);
    }

    bool push_no_duplicate(value_type const& v)
    {
        iterator i(std::upper_bound(cntnr_.begin(), cntnr_.end(), v,
                static_cast<TweakOrdering_ const&>(ord_)));
        if (i != cntnr_.begin())
        {
            if (*--i == v)
                return false;
            ++i;
        }
        cntnr_.insert(i, v);
        return true;
    }

    bool update(value_type const& v)
    {
        iterator i(std::upper_bound(cntnr_.begin(), cntnr_.end(), v,
                static_cast<TweakOrdering_ const&>(ord_)));
        if (i != cntnr_.begin())
        {
            if (*--i == v)
            {
                value_type _v(v);
                std::swap(*i, _v);
                return false;
            }
            ++i;
        }
        cntnr_.insert(i, v);
        return true;
    }


    void erase(iterator const& i)
    {
        cntnr_.erase(i);
    }

    iterator find(value_type const& v)
    {
        iterator i(std::lower_bound(begin(), end(), v,
                static_cast<TweakOrdering_ const&>(ord_)));
        return i != end() && *i == v ? i: end();
    }

    const_iterator find(value_type const& v) const
    {
        const_iterator i(std::lower_bound(begin(), end(), v,
                static_cast<TweakOrdering_ const&>(ord_)));
        return i != end() && *i == v ? i: end();
    }

    reverse_iterator rfind(value_type const& v)
    {
        reverse_iterator i(std::upper_bound(rbegin(), rend(), v,
                compose_binary(std::logical_not<bool>(),
                    static_cast<TweakOrdering_ const&>(ord_))));
        return i != rend() && *i == v ? i: rend();
    }

    const_reverse_iterator rfind(value_type const& v) const
    {
        const_reverse_iterator i(std::upper_bound(rbegin(), rend(), v,
                compose_binary(std::logical_not<bool>(),
                    static_cast<TweakOrdering_ const&>(ord_))));
        return i != rend() && *i == v ? i: rend();
    }

    size_type erase(value_type const& v)
    {
        iterator e(cntnr_.end());
        std::pair<iterator, iterator> i(std::equal_range(cntnr_.begin(), e, v,
                static_cast<TweakOrdering_ const&>(ord_)));
        const size_type retval(i.second - i.first);
        cntnr_.erase(i.first, i.second);
        return retval;
    }

    void clear()
    {
        cntnr_.clear();
    }

    sorted_list(typename boost::call_traits<TweakOrdering_>::param_type ord): ord_(ord) {}

    sorted_list(): ord_() {}

private:
    container_type cntnr_;
    TweakOrdering_ ord_;
};

#endif /* SORTED_LIST */
