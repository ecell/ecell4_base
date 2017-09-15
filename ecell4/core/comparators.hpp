#ifndef ECELL4_COMPARATORS_HPP
#define ECELL4_COMPARATORS_HPP

namespace ecell4
{

namespace utils
{

template<typename Tfirst_, typename Tsecond_>
struct pair_first_element_unary_predicator
{
    typedef std::pair<Tfirst_, Tsecond_> element_type;

    pair_first_element_unary_predicator(const Tfirst_& target)
        : target_(target)
    {
        ; // do nothing
    }

    bool operator()(const element_type& v)
    {
        return v.first == target_;
    }

protected:

    Tfirst_ target_;
};

template<typename Tfirst_, typename Tsecond_>
struct pair_second_element_unary_predicator
{
    typedef std::pair<Tfirst_, Tsecond_> element_type;

    pair_second_element_unary_predicator(const Tsecond_& target)
        : target_(target)
    {
        ; // do nothing
    }

    bool operator()(const element_type& v)
    {
        return v.second == target_;
    }

protected:

    Tsecond_ target_;
};

template<typename Tfirst_, typename Tsecond_>
struct pair_first_element_binary_predicator
    : public std::binary_function<
        std::pair<Tfirst_, Tsecond_>, std::pair<Tfirst_, Tsecond_>, bool>
{
    typedef std::pair<Tfirst_, Tsecond_> element_type;

    bool operator()(const element_type& v1, const element_type& v2)
    {
        return v1.first == v2.first;
    }
};

template<typename Tfirst_, typename Tsecond_>
struct pair_first_element_comparator
    : public std::binary_function<
        std::pair<Tfirst_, Tsecond_>, std::pair<Tfirst_, Tsecond_>, bool>
{
    typedef std::pair<Tfirst_, Tsecond_> element_type;

    inline bool operator()(const element_type& v1, const element_type& v2)
    {
        return v1.first < v2.first;
    }
};

template<typename Tfirst_, typename Tsecond_>
struct pair_second_element_comparator
    : public std::binary_function<
        std::pair<Tfirst_, Tsecond_>, std::pair<Tfirst_, Tsecond_>, bool>
{
    typedef std::pair<Tfirst_, Tsecond_> element_type;

    inline bool operator()(const element_type& v1, const element_type& v2)
    {
        return v1.second < v2.second;
    }
};

} // utils

} // ecell4

#endif /* ECELL4_COMPARATORS_HPP */
