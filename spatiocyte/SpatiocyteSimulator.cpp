#include "SpatiocyteSimulator.hpp"


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

void SpatiocyteSimulator::save_hdf5_init(std::string filename)
{
	using namespace H5;
	this->file_ = new H5File(filename, H5F_ACC_TRUNC);
	boost::scoped_ptr<Group> group(new Group(this->file_->createGroup("/SpatiocyteWorld")));
}


} // spatiocyte

} // ecell4
