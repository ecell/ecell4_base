#ifndef ECELL4_SGFRD_SHELL_VISITOR_APPLIER
#define ECELL4_SGFRD_SHELL_VISITOR_APPLIER
#include "ShellContainer.hpp"
#include "Single.hpp"
#include "Pair.hpp"
#include "Multi.hpp"
#include "Birth.hpp"

namespace ecell4
{
namespace sgfrd
{

template<typename T_shell_container>
struct shell_visitor_applier
{
    T_shell_container& container_;

    shell_visitor_applier(T_shell_container& con) : container_(con){}

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
        std::vector<ShellID> const& sids = dom.shell_ids();
        for(std::vector<ShellID>::const_iterator i(sids.begin()), e(sids.end());
                i != e; ++i)
        {
            boost::apply_visitor(f, container_.get_shell(*i));
        }
        return;
    }

    template<typename Functor>
    typename boost::enable_if<
        boost::is_same<bool, typename Functor::result_type>, bool>::type
    operator()(Functor& f, const Multi& dom)
    {
        std::vector<ShellID> const& sids = dom.shell_ids();
        for(std::vector<ShellID>::const_iterator i(sids.begin()), e(sids.end());
                i != e; ++i)
        {
            if(Functor::eval_manner::is_resolved(
                        boost::apply_visitor(f, container_.get_shell(*i))))
            {
                return Functor::eval_manner::value;
            }
        }
        return !Functor::eval_manner::value;
    }

    template<typename Functor>
    typename boost::enable_if<
        boost::is_same<void, typename Functor::result_type>, void>::type
    operator()(Functor&, const Birth&)
    {
        // Birth domain does not have any shell. do nothing.
        return ;
    }

    template<typename Functor>
    typename boost::disable_if<
        boost::is_same<void, typename Functor::result_type>, void>::type
    operator()(Functor&, const Birth&)
    {
        // Birth domain does not have any shell. do nothing.
        // if the return value is not void, return a default value.
        return typename Functor::result_type();
    }

};

} // sgfrd
} // ecell4
#endif // ECELL4_SGFRD_SHELL_VISITOR_APPLIER
