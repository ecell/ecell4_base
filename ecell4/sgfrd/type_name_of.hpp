#ifndef ECELL4_SGFRD_TYPE_NAME_OF
#define ECELL4_SGFRD_TYPE_NAME_OF
#include <boost/type_index.hpp>
#include <string>

namespace ecell4
{
namespace sgfrd
{

template<typename T>
struct type_name_of
{
    std::string operator()() const
    {
        return boost::typeindex::type_id<T>().pretty_name();
    }
};

} // sgfrd
}// ecell4
#endif// ECELL4_SGFRD_TYPE_NAME_OF
