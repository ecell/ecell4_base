#include "BDSimulator.hpp"


namespace ecell4
{

namespace bd
{

void BDSimulator::step()
{
    {
        BDPropagator propagator(*model_, *world_, rng(), dt());
        while (propagator())
        {
            ; // do nothing here
        }
    }

    set_t(t() + dt());
    ++num_steps_;
}

bool BDSimulator::step(Real const& upto)
{
    Real const t0(t()), dt0(dt());
    Real const next_time(t0 + dt0);

    if (upto <= t0)
    {
        return false;
    }

    if (upto >= next_time)
    {
        step();
        return true;
    }
    else
    {
        set_dt(upto - t0);
        step();
        set_dt(dt0);
        return false;
    }
}

void BDSimulator::save_hdf5_init(std::string filename)
{
	using namespace H5;
	file_ = new H5File(filename, H5F_ACC_TRUNC);
	boost::scoped_ptr<Group> group (new Group(file_->createGroup( "/ParticleSpace" )));
}

void BDSimulator::save(void)
{
	using namespace H5;
	if (this->file_ == NULL)
		return;

	(*world_).save_space(file_);

}

} // bd

} // ecell4
