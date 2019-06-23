#ifndef ECELL4_TYPE_NAME_OF
#define ECELL4_TYPE_NAME_OF
#include <boost/type_index.hpp>
#include <string>

namespace ecell4
{

namespace utils
{

template<typename T>
struct type_name_of
{
    std::string operator()() const
    {
        return boost::typeindex::type_id<T>().pretty_name();
    }

    static std::string value()
    {
        return boost::typeindex::type_id<T>().pretty_name();
    }
};

} // utils
} // ecell4
#endif// ECELL4_TYPE_NAME_OF
