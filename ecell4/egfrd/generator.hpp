#ifndef ECELL4_EGFRD_GENERATOR_HPP
#define ECELL4_EGFRD_GENERATOR_HPP

#include <cstddef>
#include <stdexcept>
#include <functional>

#include <boost/utility/enable_if.hpp>
#include <boost/call_traits.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/optional.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/size_type.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/mpl/or.hpp>
#include <boost/iterator/iterator_categories.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/is_readable_iterator.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_const.hpp>
#include "utils/range.hpp"

namespace ecell4
{
namespace egfrd
{

template<typename Tgen_>
bool valid(Tgen_ const& t)
{
    return true;
}


template<typename Tgen_>
std::size_t count(Tgen_ const& t)
{
    throw std::runtime_error("generation is not limited");
}

template<typename Tgen_, typename Tpred_>
bool drop_while(Tgen_& gen, Tpred_& pred)
{
    do
    {
        if (!valid(gen))
            return false;
    } while (pred(gen()));
    return true;
}

template<typename Tgen_, typename Tpred_>
bool drop_until(Tgen_& gen, Tpred_& pred)
{
    do
    {
        if (!valid(gen))
            return false;
    } while (!pred(gen()));
    return true;
}

template<typename Tgen_>
bool cue(Tgen_& gen,
    typename boost::call_traits<typename Tgen_::result_type>::param_type val)
{
    do
    {
        if (!valid(gen))
            return false;
    } while (val != gen());

    return true;
}

template<typename Tretval_>
struct abstract_generator
{
    typedef Tretval_ result_type;

    virtual ~abstract_generator() {}

    virtual Tretval_ operator()() = 0;
};


template<typename Tretval_>
struct abstract_limited_generator: public abstract_generator<Tretval_>
{
    virtual ~abstract_limited_generator() {}

    virtual std::size_t count() const
    {
        throw std::runtime_error("indetermined");
    }

    virtual bool valid() const = 0;

    virtual Tretval_ operator()() = 0;
};


template<typename Tretval_>
bool valid(abstract_limited_generator<Tretval_> const& gen)
{
    return gen.valid();
}

template<typename Tretval_>
std::size_t count(abstract_limited_generator<Tretval_> const& gen)
{
    return gen.count();
}

template<typename Trange_,
         typename Titer_ = typename boost::range_iterator<Trange_>::type,
         typename Tresult_ = typename boost::iterator_reference<Titer_>::type,
         bool Bra_ =
            boost::is_convertible<
                typename boost::iterator_category_to_traversal<
                    typename boost::BOOST_ITERATOR_CATEGORY<Titer_>::type
                    >::type,
                boost::random_access_traversal_tag>::value>
class range_generator: public abstract_limited_generator<Tresult_>
{
    template<typename Trange, typename Titer, typename Tresult, bool Bra>
    friend bool valid(range_generator<Trange, Titer, Tresult, Bra> const& gen);
    template<typename Trange, typename Titer, typename Tresult, bool Bra>
    friend std::size_t count(range_generator<Trange, Titer, Tresult, Bra> const& gen);

public:
    typedef Titer_ range_iterator;
    typedef Tresult_ result_type;

public:
    template<typename Tanother_range_>
    range_generator(Tanother_range_ const& range)
        : i_(boost::begin(range)), end_(boost::end(range)),
          count_(ecell4::egfrd::size(range)) {}

    template<typename Tanother_range_>
    range_generator(Tanother_range_& range)
        : i_(boost::begin(range)), end_(boost::end(range)),
          count_(ecell4::egfrd::size(range)) {}

    range_generator(range_iterator const& begin, range_iterator const& end)
        : i_(begin), end_(end), count_(ecell4::egfrd::size(std::make_pair(begin, end))) {}

    virtual ~range_generator() {}

    virtual result_type operator()()
    {
        --count_;
        return *i_++;
    }

    virtual std::size_t count() const
    {
        return count_;
    }

    virtual bool valid() const
    {
        return i_ != end_;
    }

private:
    range_iterator i_, end_;
    std::size_t count_;
};

template<typename Trange_, typename Titer_, typename Tresult_>
class range_generator<Trange_, Titer_, Tresult_, false>
    : public abstract_limited_generator<Tresult_>
{
    template<typename Trange, typename Titer, typename Tresult, bool Bra>
    friend bool valid(range_generator<Trange, Titer, Tresult, Bra> const& gen);
    template<typename Trange, typename Titer, typename Tresult, bool Bra>
    friend std::size_t count(range_generator<Trange, Titer, Tresult, Bra> const& gen);

public:
    typedef Titer_ range_iterator;
    typedef Tresult_ result_type;

public:
    template<typename Tanother_range_>
    range_generator(Tanother_range_ const& range)
        : i_(boost::begin(range)), end_(boost::end(range)) {}

    template<typename Tanother_range_>
    range_generator(Tanother_range_& range)
        : i_(boost::begin(range)), end_(boost::end(range)) {}

    range_generator(range_iterator const& begin, range_iterator const& end)
        : i_(begin), end_(end) {}

    template<typename Tanother_range_>
    range_generator(Tanother_range_ const& range, std::size_t count)
        : i_(boost::begin(range)), end_(boost::end(range)),
          count_(count) {}

    template<typename Tanother_range_>
    range_generator(Tanother_range_& range, std::size_t count)
        : i_(boost::begin(range)), end_(boost::end(range)),
          count_(count) {}

    range_generator(range_iterator const& begin, range_iterator const& end,
                    std::size_t count)
        : i_(begin), end_(end), count_(count) {}

    virtual ~range_generator() {}

    virtual result_type operator()()
    {
        if (count_.is_initialized())
        {
            --boost::get(count_);
        }
        return *i_++;
    }

    virtual std::size_t count() const
    {
        if (count_.is_initialized())
        {
            return boost::get(count_);
        }
        throw std::runtime_error("count not given through the constructor");
    }

    virtual bool valid() const
    {
        return i_ != end_;
    }

private:
    range_iterator i_, end_;
    boost::optional<std::size_t> count_;
};

template<bool, typename T_>
inline abstract_limited_generator<typename boost::iterator_reference<typename boost::range_iterator<T_>::type>::type >*
make_range_generator(T_& range)
{
    return new range_generator<T_, typename boost::range_iterator<T_>::type, typename boost::iterator_reference<typename boost::range_iterator<T_>::type>::type>(range);
}

template<typename Tresult_, typename T_>
inline abstract_limited_generator<Tresult_>*
make_range_generator(T_& range)
{
    return new range_generator<T_, typename boost::range_iterator<T_>::type, Tresult_>(range);
}

template<bool, typename T_>
inline abstract_limited_generator<typename boost::iterator_reference<typename boost::range_const_iterator<T_>::type>::type >*
make_range_generator(T_ const& range)
{
    return new range_generator<T_, typename boost::range_const_iterator<T_>::type, typename boost::iterator_reference<typename boost::range_const_iterator<T_>::type>::type>(range);
}

template<typename Tresult_, typename T_>
inline abstract_limited_generator<Tresult_>*
make_range_generator(T_ const& range)
{
    return new range_generator<T_, typename boost::range_const_iterator<T_>::type, Tresult_>(range);
}

} // egfrd
} // ecell4
#endif /* GENERATOR_HPP */
