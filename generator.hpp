#ifndef GENERATOR_HPP
#define GENERATOR_HPP

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
#include <boost/iterator/iterator_categories.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/is_readable_iterator.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_const.hpp>
#include "utils/reset.hpp"
#include "utils/reference_or_instance.hpp"

template<typename Tgen_>
bool valid(Tgen_ const& t)
{
    return true;
}


template<typename Tgen_>
std::size_t count(Tgen_ const& t)
{
    throw std::domain_error("generation is not limited");
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
        throw std::domain_error("indetermined");
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

template<typename Tgen_, typename Tpointer_>
class ptr_generator
{
public:
    typedef Tgen_ generator_type;
    typedef Tpointer_ pointer_type;
    typedef typename generator_type::result_type result_type;

public:
    ptr_generator(Tpointer_ const& impl): impl_(impl)
    {
        BOOST_ASSERT(&*impl_);
    }

    explicit ptr_generator(Tpointer_& impl): impl_(impl)
    {
        BOOST_ASSERT(&*impl_);
    }

    ptr_generator(ptr_generator const& that)
        : impl_(const_cast<ptr_generator&>(that).impl_) {}

    explicit ptr_generator(ptr_generator& that): impl_(that.impl_) {}

    bool valid() const
    {
        return ::valid(*impl_);
    }

    std::size_t count() const
    {
        return ::count(*impl_);
    }

    result_type operator()()
    {
        return (*impl_)();
    }

    Tpointer_ const& ptr() const
    {
        return impl_;
    }

private:
    ptr_generator const& operator=(ptr_generator const& rhs)
    {
        impl_ = rhs.impl_;
        return *this;
    }

private:
    Tpointer_ impl_;
};

template<typename Tgen, typename Tpointer>
bool valid(ptr_generator<Tgen, Tpointer> const& gen)
{
    return gen.valid();
}

template<typename Tgen, typename Tpointer>
bool count(ptr_generator<Tgen, Tpointer> const& gen)
{
    return gen.count();
}


template<typename Trange_,
         typename Titer_ = typename boost::range_iterator<Trange_>::type,
         typename Tresult_ = typename boost::iterator_reference<Titer_>::type,
         bool Bra_ = boost::is_same<typename boost::BOOST_ITERATOR_CATEGORY<Titer_>::type, std::random_access_iterator_tag>::value >
class range_generator: public abstract_limited_generator<Tresult_>
{
    template<typename Trange, typename Titer, typename Tresult, bool Bra>
    friend bool ::valid(range_generator<Trange, Titer, Tresult, Bra> const& gen);
    template<typename Trange, typename Titer, typename Tresult, bool Bra>
    friend std::size_t ::count(range_generator<Trange, Titer, Tresult, Bra> const& gen);

public:
    typedef Titer_ range_iterator;
    typedef Tresult_ result_type;

public:
    template<typename Tanother_range_>
    range_generator(Tanother_range_ const& range)
        : i_(boost::begin(range)), end_(boost::end(range)),
          count_(boost::size(range)) {}

    template<typename Tanother_range_>
    range_generator(Tanother_range_& range)
        : i_(boost::begin(range)), end_(boost::end(range)),
          count_(boost::size(range)) {}

    range_generator(range_iterator const& begin, range_iterator const& end)
        : i_(begin), end_(end), count_(boost::size(begin, end)) {}

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
    friend bool ::valid(range_generator<Trange, Titer, Tresult, Bra> const& gen);
    template<typename Trange, typename Titer, typename Tresult, bool Bra>
    friend std::size_t ::count(range_generator<Trange, Titer, Tresult, Bra> const& gen);

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
        throw std::domain_error("indetermined");
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

template<typename Tgen_, typename Tfun_, typename Tpointer_ = Tgen_*>
struct transform_generator: public abstract_limited_generator<typename Tfun_::result_type>
{
    typedef Tgen_ generator_type;
    typedef Tpointer_ pointer_type;
    typedef typename boost::remove_reference<Tfun_>::type transformer_type;
    typedef typename transformer_type::result_type result_type;

    virtual ~transform_generator() {}

    virtual std::size_t count() const
    {
        return ::count(*gen_);
    }

    virtual bool valid() const
    {
        return ::valid(*gen_);
    }

    virtual result_type operator()()
    {
        return fun_((*gen_)());
    }

    transformer_type const& functor() const
    {
        return fun_;
    }
 
    transformer_type& functor()
    {
        return fun_;
    }

    transform_generator(Tgen_* gen,
                        typename boost::call_traits<Tfun_>::param_type fun)
        : gen_(gen), fun_(fun) {}

public:
    Tfun_ fun_;
    Tpointer_ gen_;
};

template<typename Tgen_, typename Tfun_>
inline transform_generator<Tgen_, Tfun_, boost::shared_ptr<Tgen_> >
make_transform_generator(Tgen_* gen, Tfun_ const& fun)
{
    return transform_generator<Tgen_, Tfun_, boost::shared_ptr<Tgen_> >(gen, fun);
}

template<typename Tgen_, typename Tpointer_ = Tgen_*>
class generator_iterator:
    public boost::iterator_facade<
        generator_iterator<Tgen_, Tpointer_>,
        typename boost::remove_reference<typename Tgen_::result_type>::type,
        boost::single_pass_traversal_tag,
        typename Tgen_::result_type>
{
    typedef boost::iterator_facade<
        generator_iterator<Tgen_, Tpointer_>,
        typename boost::remove_reference<typename Tgen_::result_type>::type,
        boost::single_pass_traversal_tag,
        typename Tgen_::result_type> base_type;

public:
    typedef Tgen_ generator_type;
    typedef Tpointer_ pointer_type;
    typedef typename generator_type::result_type reference;

public:
    generator_iterator(): gen_(), advanced_(false) {}

    generator_iterator(Tpointer_ const& gen)
        : gen_(valid(*gen) ? gen: Tpointer_()),
          advanced_(false) {}

    void fetch()
    {
        if (gen_ && !advanced_)
        {
            if (valid(*gen_))
            {
                last_ = (*gen_)();
            }
            else
            {
                ::reset(gen_);
                last_.reset();
            }
            advanced_ = true;
        }
    }

    void increment()
    {
        fetch();
        advanced_ = false;
    }

    reference dereference() const
    {
        const_cast<generator_iterator*>(this)->fetch();
        return last_.get();
    }

    bool equal(generator_iterator const& rhs) const
    {
        const_cast<generator_iterator*>(this)->fetch();
        return (!gen_ && !rhs.gen_) || (gen_ && rhs.gen_ && *gen_ == *rhs.gen_);
    }

    typename base_type::difference_type distance_to(
            generator_iterator const& rhs) const
    {
        return (gen_ ? ::count(*gen_): 0) - (rhs.gen_ ? ::count(*rhs.gen_): 0);
    }

protected:
    Tpointer_ gen_;
    bool advanced_;
    boost::optional<reference> last_;
};

template<typename Tgen_, typename Tpointer_>
class generator_range
{
    template<typename T_> friend typename std::iterator_traits<T_>::difference_type std::distance(T_, T_);

public:
    typedef Tgen_ generator_type;
    typedef typename boost::remove_reference<typename generator_type::result_type>::type value_type;
    typedef typename generator_type::result_type reference;
    typedef const reference const_reference;

    typedef generator_iterator<Tgen_, Tpointer_> iterator;
    typedef iterator const_iterator;

public:
    iterator begin() const
    {
        return iterator(gen_);
    }

    iterator end() const
    {
        return iterator();
    }

    generator_range(Tpointer_ const& gen): gen_(gen) {}

private:
    Tpointer_ gen_;
};

namespace std {


} // namespace std

template<typename Tgen, typename Tpointer>
inline generator_range<Tgen, Tpointer>
make_generator_range(Tpointer const& gen)
{
    return generator_range<Tgen, Tpointer>(gen);
}

#endif /* GENERATOR_HPP */
