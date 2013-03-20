#include "BDSimulator.hpp"

#include <boost/scoped_array.hpp>


namespace ecell4
{

namespace bd
{

void BDSimulator::step()
{
    {
        BDPropagator propagator(*model_, *world_, *rng(), dt());
        while (propagator())
        {
            ; // do nothing here
        }
    }

    set_t(t() + dt());
    ++num_steps_;
}

bool BDSimulator::step(const Real& upto)
{
    const Real t0(t()), dt0(dt()), tnext(next_time());

    if (upto <= t0)
    {
        return false;
    }

    if (upto >= tnext)
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
	this->file_ = new H5File(filename, H5F_ACC_TRUNC);
	boost::scoped_ptr<Group> group
			(new Group(this->file_->createGroup("/BDWorld")));
}

void BDSimulator::save_hdf5(void)
{
	using namespace H5;
	const BDWorld::particle_container_type &particles = this->world_->particles();
	int length = particles.size();
	boost::scoped_array<h5_particles> h5_p(new h5_particles[length]);
	boost::scoped_array<h5_particles_index> h5_index(new h5_particles_index[length]);

	// Create Data Set.
    for (unsigned int i(0); i < length; ++i)
    {
        h5_p[i].h5_particle_id = particles[i].first;
        h5_p[i].h5_particle_position[0] = particles[i].second.position()[0];
        h5_p[i].h5_particle_position[1] = particles[i].second.position()[1];
        h5_p[i].h5_particle_position[2] = particles[i].second.position()[2];

		h5_index[i].h5_particle_id = particles[i].first;
		std::strcpy(h5_index[i].h5_particle_name, particles[i].
						second.species().name().c_str());
    }

	// Define Structure Type.
	// 	1. Positions
    CompType mtype(sizeof(h5_particles));
    mtype.insertMember(MEMBER1, HOFFSET(h5_particles, h5_particle_id),
                       PredType::NATIVE_INT);
	const hsize_t dims[] = {3};
    mtype.insertMember(MEMBER2, HOFFSET(h5_particles, h5_particle_position),
                       ArrayType(PredType::NATIVE_DOUBLE, 1, dims));
    hsize_t dim[] = {particles.size()};
    DataSpace space(1, dim);

	//	2. Tables between Id and the Name of Species.
	CompType mtype_index(sizeof(h5_particles_index));
	mtype_index.insertMember(MEMBER1, HOFFSET(h5_particles_index, h5_particle_id),
					PredType::NATIVE_INT);
	mtype_index.insertMember("Species", HOFFSET(h5_particles_index, h5_particle_name),
					StrType(PredType::C_S1, 32));


	// Create Group that represents t.
	std::ostringstream ost_hdf5path;
	boost::scoped_ptr<Group> parent_group
			(new Group (this->file_->openGroup("/BDWorld")));
	ost_hdf5path << "/BDWorld/" << this->t();
	boost::scoped_ptr<Group> group
			(new Group(parent_group->createGroup(ost_hdf5path.str())));

	// Create Dataset and Write on HDF5 File.
    DataSet *dataset = new DataSet(this->file_->createDataSet(ost_hdf5path.str() + "/Particles" , mtype, space));
    dataset->write(h5_p.get(), mtype);

    DataSet *dataset_index = new DataSet(this->file_->createDataSet(ost_hdf5path.str() + "/index" , mtype_index, space));
	dataset_index->write(h5_index.get(), mtype_index);

    delete dataset;
	delete dataset_index;

}


} // bd

} // ecell4
