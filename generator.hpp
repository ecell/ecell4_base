#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <cstddef>
#include <stdexcept>
#include <functional>
#include <boost/range/value_type.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/size_type.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/iterator_range.hpp>

template<typename Tgen_>
bool valid(Tgen_ const& t)
{
    return true;
}


template<typename Tgen_>
std::size_t count(Tgen_ const& t)
{
    throw std::domain_error();
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

    virtual bool valid() const = 0;

    virtual Tretval_ operator()() = 0;
};


template<typename Tretval_>
bool valid(abstract_limited_generator<Tretval_> const& gen)
{
    return gen.valid();
}


template<typename Trange_>
class range_generator
    : public abstract_limited_generator<typename boost::range_value<Trange_>::type>
{
    template<typename T_> friend bool ::valid(range_generator<T_> const& gen);

private:
    typedef typename boost::range_iterator<Trange_>::type range_iterator;
public:
    typedef typename boost::range_value<Trange_>::type result_type;

public:
    range_generator(Trange_ const& range)
        : i_(boost::begin(range)), end_(boost::end(range)) {}

    range_generator(std::pair<range_iterator, range_iterator> const& range)
        : i_(boost::begin(range)), end_(boost::end(range)) {}

    range_generator(boost::iterator_range<range_iterator> const& range)
        : i_(boost::begin(range)), end_(boost::end(range)) {}

    range_generator(range_iterator const& begin, range_iterator const& end)
        : i_(begin), end_(end) {}

    virtual ~range_generator() {}

    virtual result_type operator()()
    {
        return *i_++;
    }

    typename boost::range_size<Trange_>::type count()
    {
        return boost::size(std::make_pair(i_, end_));
    }

    virtual bool valid() const
    {
        return i_ != end_;
    }

private:
    range_iterator i_, end_;
};

template<typename T_>
inline abstract_limited_generator<typename T_::value_type>*
make_range_generator(T_ const& range)
{
    return new range_generator<T_>(range);
}

template<typename Trange_>
typename boost::range_size<Trange_>::type count(range_generator<Trange_> const& gen)
{
    return gen.size();
}

#endif /* GENERATOR_HPP */
