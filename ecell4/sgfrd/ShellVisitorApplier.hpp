#ifndef ECELL4_SGFRD_SHELL_VISITOR_APPLIER
#define ECELL4_SGFRD_SHELL_VISITOR_APPLIER
#include "ShellContainer.hpp"
#include "Single.hpp"
#include "Pair.hpp"
#include "Multi.hpp"

namespace ecell4
{
namespace sgfrd
{

struct minimal_eval_or : boost::true_type
{
    static bool is_resolved(const bool v) {return v;}
};

struct minimal_eval_and : boost::false_type
{
    static bool is_resolved(const bool v) {return !v;}
};

template<typename T_polygon_traits>
struct mutable_shell_visitor_applier
{
    ShellContainer<T_polygon_traits>& container_;

    mutable_shell_visitor_applier(ShellContainer<T_polygon_traits>& con)
        : container_(con)
    {}

    template<typename Functor>
    typename Functor::result_type
    operator()(Functor& f, const Single& dom)
    {
        return boost::apply_visitor(f, container_.get_shell(dom.shell_id()));
    }

    template<typename Functor>
    typename Functor::result_type
    operator()(Functor& f, const Pair& dom)
    {
        return boost::apply_visitor(f, container_.get_shell(dom.shell_id()));
    }

    template<typename Functor>
    typename boost::enable_if<
        boost::is_same<void, typename Functor::result_type>, void>::type
    operator()(Functor& f, const Multi& dom)
    {
        std::vector<ShellID> const& sids = dom.shells();
        for(typename std::vector<ShellID>::const_iterator
            iter = sids.begin(); iter != sids.end(); ++iter)
        {
            boost::apply_visitor(f, container_.get_shell(*iter));
        }
        return;
    }

    template<typename Functor>
    typename boost::enable_if<
        boost::is_same<bool, typename Functor::result_type>, bool>::type
    operator()(Functor& f, const Multi& dom)
    {
        std::vector<ShellID> const& sids = dom.shells();
        for(typename std::vector<ShellID>::const_iterator
            iter = sids.begin(); iter != sids.end(); ++iter)
        {
            if(Functor::eval_manner::is_resolved(
                        boost::apply_visitor(f, container_.get_shell(*iter))))
            {
                return Functor::eval_manner::value;
            }
        }
        return !Functor::eval_manner::value;
    }
};

} // sgfrd
} // ecell4
#endif // ECELL4_SGFRD_SHELL_VISITOR_APPLIER
