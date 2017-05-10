#ifndef EGFRD_SORTED_LIST
#define EGFRD_SORTED_LIST

#include <algorithm>
#include <functional>
#include <boost/call_traits.hpp>
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

template<typename Tcntnr_, typename TweakOrdering_ = std::less<typename boost::range_value<Tcntnr_>::type>, typename Tholder_ = Tcntnr_>
class sorted_list
{
public:
    typedef Tcntnr_ container_type;
	typedef Tholder_ holder_type;
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
        return boost::size(static_cast<container_type const&>(cntnr_));
    }

    iterator begin()
    {
        return boost::begin(static_cast<container_type&>(cntnr_));
    }

    const_iterator begin() const
    {
        return boost::begin(static_cast<container_type const&>(cntnr_));
    }

    iterator end()
    {
        return boost::end(static_cast<container_type&>(cntnr_));
    }

    const_iterator end() const
    {
        return boost::end(static_cast<container_type const&>(cntnr_));
    }

    reverse_iterator rbegin()
    {
        return boost::rbegin(static_cast<container_type&>(cntnr_));
    }

    const_reverse_iterator rbegin() const
    {
        return boost::rbegin(static_cast<container_type const&>(cntnr_));
    }

    reverse_iterator rend()
    {
        return boost::rend(static_cast<container_type&>(cntnr_));
    }

    const_reverse_iterator rend() const
    {
        return boost::end(static_cast<container_type const&>(cntnr_));
    }

    void push(value_type const& v)
    {
        iterator i(std::upper_bound(begin(), end(), v,
                static_cast<TweakOrdering_ const&>(ord_)));
        static_cast<container_type&>(cntnr_).insert(i, v);
    }

    bool push_no_duplicate(value_type const& v)
    {
        iterator i(std::upper_bound(begin(), end(), v,
                static_cast<TweakOrdering_ const&>(ord_)));
        if (i != begin())
        {
            if (*--i == v)
                return false;
            ++i;
        }
        static_cast<container_type&>(cntnr_).insert(i, v);
        return true;
    }

    bool update(value_type const& v)
    {
        iterator i(std::upper_bound(begin(), end(), v,
                static_cast<TweakOrdering_ const&>(ord_)));
        if (i != begin())
        {
            if (*--i == v)
            {
                value_type _v(v);
                std::swap(*i, _v);
                return false;
            }
            ++i;
        }
        static_cast<container_type&>(cntnr_).insert(i, v);
        return true;
    }


    void erase(iterator const& i)
    {
        static_cast<container_type&>(cntnr_).erase(i);
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
        iterator e(end());
        std::pair<iterator, iterator> i(std::equal_range(begin(), e, v,
                static_cast<TweakOrdering_ const&>(ord_)));
        const size_type retval(i.second - i.first);
        static_cast<container_type&>(cntnr_).erase(i.first, i.second);
        return retval;
    }

    void clear()
    {
        static_cast<container_type&>(cntnr_).clear();
    }

	holder_type& container()
	{
		return cntnr_;
	}

	holder_type const& container() const
	{
		return cntnr_;
	}

    sorted_list(typename boost::call_traits<TweakOrdering_>::param_type ord,
			    typename boost::call_traits<holder_type>::param_type holder)
		: ord_(ord), cntnr_(holder) {}

    explicit sorted_list(typename boost::call_traits<holder_type>::param_type holder)
		: ord_(), cntnr_(holder) {}

    explicit sorted_list(typename boost::call_traits<TweakOrdering_>::param_type ord): ord_(ord) {}

    sorted_list(): ord_() {}

private:
    TweakOrdering_ ord_;
    holder_type cntnr_;
};

#endif /* EGFRD_SORTED_LIST */
