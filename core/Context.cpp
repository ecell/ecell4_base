#include "Context.hpp"


namespace ecell4
{

bool uspmatch(const UnitSpecies& pttrn, const Species& sp)
{
    for (Species::container_type::const_iterator i(sp.begin());
        i != sp.end(); ++i)
    {
        const UnitSpecies& usp(*i);
        bool succeeded(true);
        for (UnitSpecies::container_type::const_iterator j(pttrn.begin());
            j != pttrn.end(); ++j)
        {
            if (usp.has_site((*j).first))
            {
                const UnitSpecies::site_type& site(usp.get_site((*j).first));
                if ((*j).second.first != "" && (*j).second.first != site.first)
                {
                    succeeded = false;
                    break;
                }
                else if ((*j).second.second == "" && site.second != "")
                {
                    succeeded = false;
                    break;
                }

                ; //XXX: pass
            }
            else
            {
                succeeded = false;
                break;
            }
        }

        if (succeeded)
        {
            return true;
        }
    }
    return false;
}

bool spmatch(const Species& pttrn, const Species& sp)
{
    // Context ctx;

    for (Species::container_type::const_iterator i(pttrn.begin());
        i != pttrn.end(); ++i)
    {
        const UnitSpecies& usp(*i);
        if (!uspmatch(usp, sp))
        {
            return false;
        }
    }

    return true;
}

} // ecell4
