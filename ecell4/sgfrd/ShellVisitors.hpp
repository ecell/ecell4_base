#ifndef ECELL4_SGFRD_SHELL_VISITORS
#define ECELL4_SGFRD_SHELL_VISITORS
#include "Shell.hpp"
#include "SGFRDEvent.hpp"

namespace ecell4
{
namespace sgfrd
{

struct domain_id_setter : boost::static_visitor<void>
{
    DomainID did_;
    domain_id_setter(DomainID did) : did_(did){}

    template<typename shapeT, typename stridT>
    void operator()(Shell<shapeT, stridT>& shell) const
    {
        shell.domain_id() = did_;
    }
};

struct domain_id_getter : boost::static_visitor<DomainID>
{
    template<typename shapeT, typename stridT>
    DomainID operator()(const Shell<shapeT, stridT>& shell) const
    {
        return shell.domain_id();
    }
};






} // sgfrd
} // ecell4
#endif// ECELL4_SGFRD_SHELL_VISITORS
