#ifndef GENERATOR_HPP
#define GENERATOR_HPP

#include <cstddef>
#include <stdexcept>
#include <functional>

#include <boost/utility/enable_if.hpp>
#include <boost/call_traits.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/size.hpp>
#include <boost/range/size_type.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include <boost/iterator/is_readable_iterator.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_const.hpp>
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

template<typename Timpl_, typename Tholder_ = boost::shared_ptr<Timpl_> >
class ptr_generator
{
public:
    typedef typename Timpl_::result_type result_type;

private:
    typedef Tholder_ holder_type;

public:
    ptr_generator(Timpl_* impl): impl_(impl)
    {
        BOOST_ASSERT(impl);
    }

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

    holder_type ptr()
    {
        return impl_;
    }

private:
    ptr_generator& operator=(ptr_generator const&) { return *this; }

private:
    holder_type impl_;
    int refcount_;
};

template<typename Timpl_>
bool valid(ptr_generator<Timpl_> const& gen)
{
    return gen.valid();
}

template<typename Timpl_>
bool count(ptr_generator<Timpl_> const& gen)
{
    return gen.count();
}


template<typename Trange_, typename Titer_ = typename boost::range_iterator<Trange_>::type, typename Tresult_ = typename boost::iterator_reference<Titer_>::type>
class range_generator: public abstract_limited_generator<Tresult_>
{
    template<typename T_> friend bool ::valid(range_generator<T_> const& gen);

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
    typedef Tfun_ transformer_type;
    typedef typename Tfun_::result_type result_type;

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

    transform_generator(Tgen_* gen, Tfun_ const& fun, int)
        : gen_(gen), fun_(fun, 1) {}

    transform_generator(Tgen_* gen, Tfun_ fun)
        : gen_(gen), fun_(fun, 1) {}

    transform_generator(Tgen_* gen, Tfun_& fun)
        : gen_(gen), fun_(fun) {}

public:
    reference_or_instance<Tfun_> fun_;
    Tpointer_ gen_;
};

template<typename Tgen_, typename Tfun_>
inline transform_generator<Tgen_, Tfun_, boost::shared_ptr<Tgen_> >
make_transform_generator(Tgen_* gen, Tfun_ const& fun)
{
    return transform_generator<Tgen_, Tfun_, boost::shared_ptr<Tgen_> >(gen, fun);
}


template<typename Tgen_, typename Tpointer_ = Tgen_*>
class generator_range
{
    template<typename T_> friend typename std::iterator_traits<T_>::difference_type std::distance(T_, T_);

public:
    class iterator
    {
    public:
        typedef std::forward_iterator_tag iterator_category;
        typedef std::ptrdiff_t difference_type;
        typedef typename Tgen_::result_type reference;
        typedef typename boost::remove_reference<reference>::type value_type;
        typedef value_type* pointer;

    public:
        iterator(): gen_(0), advanced_(false) {}

        iterator(Tgen_* gen): gen_(valid(*gen) ? gen: 0), advanced_(false) {}

        iterator& operator++(int)
        {
            if (gen_ && !advanced_)
            {
                if (valid(*gen_))
                    last_ = (*gen_)();
                else
                    gen_ = 0;
            }
            advanced_ = false;
            return *this;
        }

        reference operator*() const
        {
            if (gen_ && !advanced_)
            {
                if (valid(*gen_))
                    last_ = (*gen_)();
                else
                    gen_ = 0;
            }
            advanced_ = true;
            return last_;
        }

        bool operator==(iterator const& rhs) const
        {
            return gen_ == rhs.gen_;
        }

        bool operator!=(iterator const& rhs) const
        {
            return operator!=(rhs);
        }

    protected:
        Tgen_* gen_;
        bool advanced_;
        typename boost::remove_const<value_type>::type last_;
    };

    typedef iterator const_iterator;

public:
    iterator begin() const
    {
        return iterator(&*gen_);
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

template<typename Tgen_, typename Tpointer_>
inline typename generator_range<Tgen_, Tpointer_>::difference_type
distance(generator_range<Tgen_, Tpointer_> const& lhs, generator_range<Tgen_, Tpointer_> const& rhs)
{
    return ::count(*lhs.gen_);
}

} // namespace std

template<typename Tgen_>
inline generator_range<Tgen_, boost::shared_ptr<Tgen_> >
make_generator_range(Tgen_* gen)
{
    return generator_range<Tgen_, boost::shared_ptr<Tgen_> >(
            boost::shared_ptr<Tgen_>(gen));
}

#endif /* GENERATOR_HPP */
