#include "SpatiocyteSimulator.hpp"

#include <boost/scoped_array.hpp>

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

const H5std_string SPATIOCYTE_MEMBER1("lattice_id");
//const H5std_string SPATIOCYTE_MEMBER2("positions");
const H5std_string SPATIOCYTE_MEMBER2("species_id");


namespace ecell4
{

namespace spatiocyte
{

void SpatiocyteSimulator::step()
{
    (*world_).step();

    ++num_steps_;
}

bool SpatiocyteSimulator::step(const Real& upto)
{
    const bool retval((*world_).step(upto));
    ++num_steps_;
    return retval;
}

} // spatiocyte

} // ecell4
