#ifndef ECELL4_SGFRD_VISIT_SHELL
#define ECELL4_SGFRD_VISIT_SHELL
#include <ecell4/sgfrd/Single.hpp>
#include <ecell4/sgfrd/Pair.hpp>
#include <ecell4/sgfrd/Multi.hpp>
#include <boost/utility/enable_if.hpp>

namespace ecell4
{
namespace sgfrd
{

template<typename Functor>
inline typename Functor::result_type
visit_shell(Functor& f, const Single& s)
{
    return boost::apply_visitor(f, s.shell());
}

template<typename Functor>
inline typename Functor::result_type
visit_shell(Functor& f, Single& s)
{
    return boost::apply_visitor(f, s.shell());
}

template<typename Functor>
inline typename Functor::result_type
visit_shell(Functor& f, const Pair& s)
{
    return boost::apply_visitor(f, s.shell());
}

template<typename Functor>
inline typename Functor::result_type
visit_shell(Functor& f, Pair& s)
{
    return boost::apply_visitor(f, s.shell());
}

template<typename Functor>
inline typename Functor::result_type
visit_shell(Functor& f, const Multi& s, std::size_t idx)
{
    return boost::apply_visitor(f, s.shells().at(idx));
}

template<typename Functor>
inline typename Functor::result_type
visit_shell(Functor& f, Multi& s, std::size_t idx)
{
    return boost::apply_visitor(f, s.shells().at(idx));
}

template<typename Functor>
typename boost::enable_if<
    boost::is_same<typename Functor::result_type, void>, void>::type
visit_shells(Functor& f, const Multi& s)
{
    for(typename Multi::shell_container_type::const_iterator
        iter(s.shells().begin()), end(s.shells().end()); iter != end; ++iter)
    {
        boost::apply_visitor(f, *iter);
    }
    return;
}

template<typename Functor>
typename boost::enable_if<
    boost::is_same<typename Functor::result_type, void>, void>::type
visit_shells(Functor& f, Multi& s)
{
    for(typename Multi::shell_container_type::iterator
        iter(s.shells().begin()), end(s.shells().end()); iter != end; ++iter)
    {
        boost::apply_visitor(f, *iter);
    }
    return;
}

} // sgfrd
} // ecell4
#endif// ECELL4_SGFRD_VISIT_SHELL
