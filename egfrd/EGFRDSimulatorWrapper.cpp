#include "EGFRDSimulatorWrapper.hpp"

#include <boost/scoped_array.hpp>

#define MEMBER1 "particle_id_lot"
#define MEMBER2 "particle_id_serial"
#define MEMBER3 "positions"


namespace ecell4
{

namespace egfrd
{

void EGFRDSimulatorWrapper::step()
{
    if ((*world_).num_particles() == 0)
    {
        ; // increment time
        return;
    }

    (*sim_).step();
    (*world_).set_t((*sim_).t());
}

bool EGFRDSimulatorWrapper::step(const Real& upto)
{
    if ((*world_).num_particles() == 0)
    {
        ; // increment time
        return true; // this should be false
    }

    const bool retval((*sim_).step(upto));
    (*world_).set_t((*sim_).t());
    return retval;
}

void EGFRDSimulatorWrapper::save_hdf5_init(std::string filename)
{
    using namespace H5;
    this->file_ = new H5File(filename, H5F_ACC_TRUNC);
    boost::scoped_ptr<Group> group(
        new Group(this->file_->createGroup("/EGFRDWorld")));
}

void EGFRDSimulatorWrapper::save_hdf5(void)
{
    using namespace H5;
    typedef std::vector<std::pair<ParticleID, Particle> > particle_container_type;
    const particle_container_type &particles = this->world_->list_particles();

    unsigned int length = particles.size();
    boost::scoped_array<h5_particles> h5_p(new h5_particles[length]);
    boost::scoped_array<h5_particles_index>
        h5_index(new h5_particles_index[length]);

    for (unsigned int i(0); i < length; ++i)
    {
        h5_p[i].h5_particle_id_lot = particles[i].first.lot();
        h5_p[i].h5_particle_id_serial = particles[i].first.serial();
        h5_p[i].h5_particle_position[0] = particles[i].second.position()[0];
        h5_p[i].h5_particle_position[1] = particles[i].second.position()[1];
        h5_p[i].h5_particle_position[2] = particles[i].second.position()[2];

        h5_index[i].h5_particle_id_serial = particles[i].first.serial();
        std::strcpy(h5_index[i].h5_particle_name,
                    particles[i].second.species().name().c_str());

        // h5_index[i].h5_particle_radius = particles[i].second.radius();
        // h5_index[i].h5_particle_D = particles[i].second.D();

        // h5_index[i].h5_particle_radius =
        //     particles[i].second.species().get_attribute("radius");
        // h5_index[i].h5_particle_D =
        //     particles[i].second.species().get_attribute("D");

        // h5_index[i].h5_particle_radius =
        //     (this->world_->get_molecule_info(
        //          particles[i].second.species())).radius;
        // h5_index[i].h5_particle_D =
        //     (this->world_->get_molecule_info
        //          (particles[i].second.species())).D;

        // h5_index[i].h5_particle_radius = 0.0;
        // h5_index[i].h5_particle_D = 0.0;
    }

    // Define Structure Type.
    // 1. Positions
    CompType mtype(sizeof(h5_particles));
    mtype.insertMember(MEMBER1, HOFFSET(h5_particles, h5_particle_id_lot),
                       PredType::NATIVE_INT);
    mtype.insertMember(MEMBER2, HOFFSET(h5_particles, h5_particle_id_serial),
                       PredType::NATIVE_INT);
    const hsize_t dims[] = {3};
    mtype.insertMember(MEMBER3, HOFFSET(h5_particles, h5_particle_position),
                       ArrayType(PredType::NATIVE_DOUBLE, 1, dims));
    hsize_t dim[] = {particles.size()};
    DataSpace space(1, dim);

    // 2. Tables between Id and the Name of Species.
    CompType mtype_index(sizeof(h5_particles_index));
    mtype_index.insertMember(
        MEMBER2, HOFFSET(h5_particles_index, h5_particle_id_serial),
        PredType::NATIVE_INT);
    mtype_index.insertMember(
        "Species", HOFFSET(h5_particles_index, h5_particle_name),
        StrType(PredType::C_S1, 32));
    mtype_index.insertMember(
        "Radius", HOFFSET(h5_particles_index, h5_particle_radius),
        PredType::NATIVE_DOUBLE);
    mtype_index.insertMember(
        "D", HOFFSET(h5_particles_index, h5_particle_D),
        PredType::NATIVE_DOUBLE);

    // Create Group that represents t.
    std::ostringstream ost_hdf5path;
    boost::scoped_ptr<Group> parent_group(
        new Group (this->file_->openGroup("/EGFRDWorld")));
    ost_hdf5path << "/EGFRDWorld/" << this->t();
    boost::scoped_ptr<Group> group(
        new Group(parent_group->createGroup(ost_hdf5path.str())));

    // Set Attribute
    const double t_value = this->t();
    FloatType doubleType(PredType::IEEE_F64LE);

    // Create Dataset and Write on HDF5 File.
    DataSet *dataset = new DataSet(
        this->file_->createDataSet(
            ost_hdf5path.str() + "/Particles" , mtype, space));
    Attribute attr = dataset->createAttribute(
        "t", doubleType, DataSpace(H5S_SCALAR));
    attr.write(doubleType, &t_value);
    dataset->write(h5_p.get(), mtype);

    DataSet *dataset_index = new DataSet(
        this->file_->createDataSet(
            ost_hdf5path.str() + "/index" , mtype_index, space));
    Attribute attr_index = dataset_index->createAttribute
        ("t", doubleType, DataSpace(H5S_SCALAR));
    attr_index.write(doubleType, &t_value);
    dataset_index->write(h5_index.get(), mtype_index);

    delete dataset;
    delete dataset_index;
}

} // egfrd

} // ecell4
