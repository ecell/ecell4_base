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

void SpatiocyteSimulator::save_hdf5_init(std::string filename)
{
    using namespace H5;
    this->file_ = new H5File(filename, H5F_ACC_TRUNC);
    boost::scoped_ptr<Group> group(
        new Group(this->file_->createGroup("/SpatiocyteWorld")));
}

void SpatiocyteSimulator::save_hdf5(void)
{
    using namespace H5;
    const SpatiocyteWorld::species_container_type& species =
        this->world_->species();
    unsigned int species_size = species.size();
    int lattice_size = 0;

    for (unsigned int i(0); i < species_size; ++i)
    {
        lattice_size += this->world_->list_particles(species[i].first).size();
    }

    boost::scoped_array<h5_lattice> h5_lattice_p(new h5_lattice[lattice_size]);

    // Create Data Set.

    for (unsigned int i(0); i < species_size; ++i)
    {
        int tmp = 0;
        for (unsigned int j(0);
             j < this->world_->list_particles(species[i].first).size(); ++j)
        {
            h5_lattice_p[tmp + j].h5_lattice_id =
                static_cast<int>(
                    this->world_->list_particles(
                        species[i].first)[j].first.serial());

            std::strcpy(h5_lattice_p[tmp + i].h5_species_id,
                        species[i].first.c_str());
        }
        tmp += this->world_->list_particles(species[i].first).size();
    }

    // Define Structure Type
    CompType mtype(sizeof(h5_lattice));
    mtype.insertMember(SPATIOCYTE_MEMBER1,
                       HOFFSET(h5_lattice, h5_lattice_id),
                       PredType::NATIVE_INT);
    // const hsize_t dims[] = {3};
    // mtype.insertMember(SPATIOCYTE_MEMBER2,
    //                    HOFFSET(h5_particles, h5_particle_position),
    //                    ArrayType(PredType::NATIVE_DOUBLE, 1, dims));
    mtype.insertMember(SPATIOCYTE_MEMBER2,
                       HOFFSET(h5_lattice, h5_species_id),
                       StrType(PredType::C_S1, 32));
    hsize_t dim[] = {lattice_size};
    DataSpace space(1, dim);

    // Create group for each t
    std::ostringstream ost_hdf5path;
    boost::scoped_ptr<Group> parent_group(
        new Group(this->file_->openGroup("/SpatiocyteWorld")));
    ost_hdf5path << "/SpatiocyteWorld/" << this->t();
    boost::scoped_ptr<Group> group(
        new Group(parent_group->createGroup(ost_hdf5path.str())));

    // Create Dataset and Write it to HDF5
    H5::DataSet *dataset = new DataSet(
        this->file_->createDataSet(
            ost_hdf5path.str() + "/Particles", mtype, space));
    dataset->write(h5_lattice_p.get(), mtype);

    delete dataset;
}

} // spatiocyte

} // ecell4
