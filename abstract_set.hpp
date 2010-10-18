#ifndef ABSTRACT_SET_HPP
#define ABSTRACT_SET_HPP

#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/value_type.hpp>
#include <algorithm>

template<typename T_>
struct collection_value: public boost::range_value<T_>
{
};

template<typename T_>
struct inserter: public std::unary_function<typename collection_value<T_>::type, bool>
{
    typedef T_ set_type;
    typedef typename collection_value<set_type>::type argument_type;

    inserter(set_type& set): set_(set) {}

    bool operator()(argument_type const& v)
    {
        set_.push_back(v); 
        return true;
    }

private:
    set_type& set_;
};

template<typename Tval_, typename Tcompare_, typename Talloc_>
struct inserter<std::set<Tval_, Tcompare_, Talloc_> >: public std::unary_function<typename collection_value<std::set<Tval_, Tcompare_, Talloc_> >::type, bool>
{
    typedef std::set<Tval_, Tcompare_, Talloc_> set_type;
    typedef typename collection_value<set_type>::type argument_type;

    inserter(set_type& set): set_(set) {}

    bool operator()(argument_type const& v)
    {
        return set_.insert(v).second; 
    }

private:
    set_type& set_;
};

template<typename Tkey_, typename Tval_, typename Tcompare_, typename Talloc_>
struct inserter<std::map<Tkey_, Tval_, Tcompare_, Talloc_> >: public std::unary_function<typename collection_value<std::map<Tkey_, Tval_, Tcompare_, Talloc_> >::type, bool>
{
    typedef std::map<Tkey_, Tval_, Tcompare_, Talloc_> set_type;
    typedef typename collection_value<set_type>::type argument_type;

    inserter(set_type& set): set_(set) {}

    bool operator()(argument_type const& v)
    {
        return set_.insert(v).second; 
    }

private:
    set_type& set_;
};

template<typename T_>
inline bool contains(T_ const& s, typename collection_value<T_>::type const& v)
{
    typename boost::range_const_iterator<T_>::type e(boost::end(s));
    return e != std::find(boost::begin(s), e, v);
}

template<typename T_>
inline bool insert(T_& s, typename collection_value<T_>::type const& v)
{
    return inserter<T_>(s)(v);
}

template<typename T1, typename T2, typename Tr>
inline void difference(T1 const& r1, T2 const& r2, Tr const& result)
{
    std::set_difference(
        boost::begin(r1), boost::end(r1),
        boost::begin(r2), boost::end(r2),
        result);
}

#endif /* ABSTRACT_SET_HPP */
