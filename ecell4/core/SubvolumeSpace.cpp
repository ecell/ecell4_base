#include "SubvolumeSpace.hpp"
#include "Context.hpp"

namespace ecell4
{

Integer SubvolumeSpaceVectorImpl::num_molecules(
    const Species& sp, const coordinate_type& c) const
{
    SpeciesExpressionMatcher sexp(sp);
    const cell_type& cell(matrix_[c]);
    Integer retval(0);
    for (cell_type::const_iterator i(cell.begin());
        i != cell.end(); ++i)
    {
        if (sexp.match((*i).first))
        {
            do
            {
                retval += (*i).second;
            } while (sexp.next());
        }
    }
    return retval;
}

Integer SubvolumeSpaceVectorImpl::num_molecules_exact(
    const Species& sp, const coordinate_type& c) const
{
    const cell_type& cell(matrix_[c]);
    cell_type::const_iterator i(cell.find(sp));
    if (i == cell.end())
    {
        return 0;
    }
    return (*i).second;
}

void SubvolumeSpaceVectorImpl::reserve_species(
    const Species& sp, const coordinate_type& c)
{
    cell_type& cell(matrix_[c]);
    cell_type::const_iterator i(cell.find(sp));
    if (i != cell.end())
    {
        throw AlreadyExists("Species already exists");
    }
    cell.insert(std::make_pair(sp, 0));
    species_.push_back(sp);
}

void SubvolumeSpaceVectorImpl::add_molecules(
    const Species& sp, const Integer& num, const coordinate_type& c)
{
    cell_type& cell(matrix_[c]);
    cell_type::iterator i(cell.find(sp));
    if (i == cell.end())
    {
        reserve_species(sp, c);
        i = cell.find(sp);
    }
    (*i).second += num;
}

void SubvolumeSpaceVectorImpl::remove_molecules(
    const Species& sp, const Integer& num, const coordinate_type& c)
{
    cell_type& cell(matrix_[c]);
    cell_type::iterator i(cell.find(sp));
    if (i == cell.end())
    {
        std::ostringstream message;
        message << "Speices [" << sp.serial() << "] not found";
        throw NotFound(message.str());
    }

    if ((*i).second < num)
    {
        std::ostringstream message;
        message << "The number of molecules cannot be negative. [" << sp.serial() << "]";
        throw std::invalid_argument(message.str());
    }

    (*i).second -= num;
}

} // ecell4
