#ifndef UTILS_BASE_TYPE_WALKER_HPP
#define UTILS_BASE_TYPE_WALKER_HPP

#include <boost/call_traits.hpp>

template<typename Twalker_, typename T_, typename TwalkerHolder_ = Twalker_ const&>
class base_type_walker
{
public:
    typedef Twalker_ walker_type;
    typedef TwalkerHolder_ walker_holder_type;

    base_type_walker(typename boost::call_traits<walker_holder_type>::param_type walker): walker(walker) {}

    static const bool has_base_class = sizeof(typename T_::base_type**) == sizeof(void*);

    void operator()()
    {
        return this->operator()<T_>(0);
    }

private:
    template<typename T>
    void operator()(typename T::base_type*) const
    {
        walker.template operator()<T>();
        base_type_walker<Twalker_, typename T::base_type>(walker).operator()();
    }

    template<typename T>
    void operator()(...) const
    {
        walker.template operator()<T>();
    }

public:
    walker_holder_type walker;
};

#endif /* UTILS_BASE_TYPE_WALKER_HPP */
