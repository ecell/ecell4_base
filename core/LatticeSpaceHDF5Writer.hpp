#ifndef __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP
#define __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP

#include <cstring>
#include <sstream>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/lexical_cast.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"
#include "Species.hpp"


namespace ecell4
{

struct LatticeSpaceHDF5Traits
{
    typedef struct h5_species_struct {
        uint32_t id;
        char serial[32]; // species' serial may exceed the limit
    } h5_species_struct;

    typedef struct h5_voxel_struct {
        int lot;
        int serial;
        uint32_t sid;
        int64_t coord;
        double D;
    } h5_voxel_struct;

    static H5::CompType get_voxel_comp_type()
    {
        H5::CompType h5_voxel_comp_type(sizeof(h5_voxel_struct));
        h5_voxel_comp_type.insertMember(
            std::string("lot"), HOFFSET(h5_voxel_struct, lot),
            H5::PredType::NATIVE_INT);
        h5_voxel_comp_type.insertMember(
            std::string("serial"), HOFFSET(h5_voxel_struct, serial),
            H5::PredType::NATIVE_INT);
        h5_voxel_comp_type.insertMember(
            std::string("sid"), HOFFSET(h5_voxel_struct, sid),
            H5::PredType::STD_I32LE);
        h5_voxel_comp_type.insertMember(
            std::string("coordinate"), HOFFSET(h5_voxel_struct, coord),
            H5::PredType::STD_I64LE);
        h5_voxel_comp_type.insertMember(
            std::string("D"), HOFFSET(h5_voxel_struct, D),
            H5::PredType::NATIVE_DOUBLE);
        return h5_voxel_comp_type;
    }

    static H5::CompType get_species_comp_type()
    {
        H5::CompType h5_species_comp_type(sizeof(h5_species_struct));
        h5_species_comp_type.insertMember(
            std::string("id"), HOFFSET(h5_species_struct, id),
            H5::PredType::STD_I32LE);
        h5_species_comp_type.insertMember(
            std::string("serial"), HOFFSET(h5_species_struct, serial),
            H5::StrType(H5::PredType::C_S1, 32));
        return h5_species_comp_type;
    }
};

template<typename Tspace_>
void save_lattice_space(const Tspace_& space, H5::Group* root)
{
    typedef LatticeSpaceHDF5Traits traits_type;
    typedef typename traits_type::h5_species_struct h5_species_struct;
    typedef typename traits_type::h5_voxel_struct h5_voxel_struct;

    typedef std::vector<std::pair<ParticleID, Voxel> >
        voxel_container_type;

    const unsigned int num_voxels(space.num_voxels());
    const std::vector<Species> species(space.list_species());
    boost::scoped_array<h5_voxel_struct>
        h5_voxel_table(new h5_voxel_struct[num_voxels]);
    boost::scoped_array<h5_species_struct>
        h5_species_table(new h5_species_struct[species.size()]);
    unsigned int vidx(0);
    for (unsigned int i(0); i < species.size(); ++i)
    {
        const unsigned int sid(i + 1);
        h5_species_table[i].id = sid;
        std::strcpy(h5_species_table[i].serial, species[i].serial().c_str());

        const voxel_container_type& voxels(space.list_voxels(species[i]));
        for (unsigned int j(0); j < voxels.size(); ++j)
        {
            h5_voxel_table[vidx].lot = voxels[j].first.lot();
            h5_voxel_table[vidx].serial = voxels[j].first.serial();
            h5_voxel_table[vidx].sid = sid;
            h5_voxel_table[vidx].coord = voxels[j].second.coordinate();
            h5_voxel_table[vidx].D = voxels[j].second.D();
            ++vidx;
        }
    }

    const int RANK = 1;
    hsize_t dim1[] = {num_voxels};
    H5::DataSpace dataspace1(RANK, dim1);
    boost::scoped_ptr<H5::DataSet> dataset1(new H5::DataSet(
        root->createDataSet(
            "voxels", traits_type::get_voxel_comp_type(), dataspace1)));

    hsize_t dim2[] = {species.size()};
    H5::DataSpace dataspace2(RANK, dim2);
    boost::scoped_ptr<H5::DataSet> dataset2(new H5::DataSet(
        root->createDataSet(
            "species", traits_type::get_species_comp_type(), dataspace2)));

    dataset1->write(h5_voxel_table.get(), dataset1->getDataType());
    dataset2->write(h5_species_table.get(), dataset2->getDataType());

    const uint32_t space_type = static_cast<uint32_t>(Space::LATTICE);
    H5::Attribute attr_space_type(
        root->createAttribute(
            "type", H5::PredType::STD_I32LE, H5::DataSpace(H5S_SCALAR)));
    attr_space_type.write(H5::PredType::STD_I32LE, &space_type);

    const double t = space.t();
    H5::Attribute attr_t(
        root->createAttribute(
            "t", H5::PredType::IEEE_F64LE, H5::DataSpace(H5S_SCALAR)));
    attr_t.write(H5::PredType::IEEE_F64LE, &t);

    const double voxel_radius = space.voxel_radius();
    H5::Attribute attr_voxel_radius(
        root->createAttribute(
            "voxel_radius", H5::PredType::IEEE_F64LE, H5::DataSpace(H5S_SCALAR)));
    attr_voxel_radius.write(H5::PredType::IEEE_F64LE, &voxel_radius);

    const Position3 edge_lengths = space.edge_lengths();
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    H5::Attribute attr_lengths(
        root->createAttribute(
            "edge_lengths", lengths_type, H5::DataSpace(H5S_SCALAR)));
    double lengths[] = {edge_lengths[0], edge_lengths[1], edge_lengths[2]};
    attr_lengths.write(lengths_type, lengths);
}

// template<typename Tspace_>
// class LatticeSpaceHDF5Writer
// {
// public:
// 
//     typedef Tspace_ space_type;
// 
// protected:
// 
//     struct h5_species_struct
//     {
//         uint32_t id;
//         char serial[32]; // species' serial may exceed the limit
//         double radius;
//         double D;
//     };
// 
//     struct h5_particle_struct
//     {
//         uint32_t id;
//         uint32_t spid;
//         double pos[3];
//     };
// 
//     struct h5_voxel_struct
//     {
//         uint32_t sid;
//     };
// 
// public:
// 
//     LatticeSpaceHDF5Writer(const space_type& space)
//         : space_(space)
//     {
//         ;
//     }
// 
//     virtual ~LatticeSpaceHDF5Writer()
//     {
//         ;
//     }
// 
//     void save(H5::H5File* fout, const std::string& hdf5path) const
//     {
//         using namespace H5;
// 
//         std::vector<std::pair<ParticleID, Particle> >
//             particles(space_.list_particles());
//         std::vector<Species> species(space_.list_species());
// 
//         boost::scoped_array<h5_species_struct>
//             h5_species_table(new h5_species_struct[species.size()]);
//         for (unsigned int i(0); i < species.size(); ++i)
//         {
//             Species sp(species.at(i));
//             h5_species_table[i].id = i + 1;
//             std::strncpy(h5_species_table[i].serial,
//                          sp.serial().c_str(),
//                          sizeof(h5_species_table[i].serial));
//             h5_species_table[i].radius = boost::lexical_cast<double>(
//                     sp.get_attribute("radius"));
//             h5_species_table[i].D = boost::lexical_cast<double>(
//                     sp.get_attribute("D"));
//         }
// 
//         boost::scoped_array<h5_particle_struct>
//             h5_particle_table(new h5_particle_struct[particles.size()]);
//         for (unsigned int i(0); i < particles.size(); ++i)
//         {
//             const ParticleID pid(particles.at(i).first);
//             const Particle p(particles.at(i).second);
//             h5_particle_table[i].id = pid.serial(); // TODO
//             uint32_t spid(0);
//             for (unsigned int j(0); j < species.size(); ++j)
//             {
//                 if (p.species() == species.at(j))
//                 {
//                     spid = j + 1;
//                     break;
//                 }
//             }
//             h5_particle_table[i].spid = spid;
//             h5_particle_table[i].pos[0] = p.position()[0];
//             h5_particle_table[i].pos[1] = p.position()[1];
//             h5_particle_table[i].pos[2] = p.position()[2];
//         }
// 
//         hsize_t dim0[] = {3};
//         ArrayType particle_pos_type(
//             PredType::NATIVE_DOUBLE, 1, dim0);
// 
//         CompType h5_particle_comp_type(sizeof(h5_particle_struct));
//         h5_particle_comp_type.insertMember(
//             std::string("id"), HOFFSET(h5_particle_struct, id),
//             PredType::STD_I32LE);
//         h5_particle_comp_type.insertMember(
//             std::string("spid"), HOFFSET(h5_particle_struct, spid),
//             PredType::STD_I32LE);
//         h5_particle_comp_type.insertMember(
//             std::string("pos"), HOFFSET(h5_particle_struct, pos),
//             particle_pos_type);
// 
//         CompType h5_species_comp_type(sizeof(h5_species_struct));
//         h5_species_comp_type.insertMember(
//             std::string("id"), HOFFSET(h5_species_struct, id),
//             PredType::STD_I32LE);
//         h5_species_comp_type.insertMember(
//             std::string("serial"), HOFFSET(h5_species_struct, serial),
//             StrType(PredType::C_S1, 32));
//         h5_species_comp_type.insertMember(
//             std::string("radius"), HOFFSET(h5_species_struct, radius),
//             PredType::NATIVE_DOUBLE);
//         h5_species_comp_type.insertMember(
//             std::string("D"), HOFFSET(h5_species_struct, D),
//             PredType::NATIVE_DOUBLE);
// 
//         const int RANK = 1;
// 
//         // Create the group "/data"
//         std::ostringstream oss;
//         oss << hdf5path << "/data";
//         boost::scoped_ptr<Group> data_group(
//             new Group(fout->createGroup(oss.str())));
// 
//         Attribute attr_world_size(
//                 data_group->createAttribute(
//                     "world_size", PredType::IEEE_F64LE, DataSpace(H5S_SCALAR)));
//         attr_world_size.write(PredType::IEEE_F64LE, &space_.edge_lengths()[0]);
// 
//         oss << "/" << space_.t();
//         boost::scoped_ptr<Group> group(
//             new Group(data_group->createGroup(oss.str())));
// 
//         hsize_t dim1[] = {particles.size()};
//         DataSpace dataspace1(RANK, dim1);
//         boost::scoped_ptr<DataSet> dataset1(
//             new DataSet(group->createDataSet(
//                             "particles", h5_particle_comp_type,
//                             dataspace1)));
// 
//         hsize_t dim2[] = {species.size()};
//         DataSpace dataspace2(RANK, dim2);
//         boost::scoped_ptr<DataSet> dataset2(
//             new DataSet(fout->createDataSet(
//                             hdf5path + "/species" , h5_species_comp_type,
//                             dataspace2)));
// 
//         dataset1->write(h5_particle_table.get(), h5_particle_comp_type);
//         dataset2->write(h5_species_table.get(), h5_species_comp_type);
// 
//         const double t = space_.t();
//         Attribute attr_t(
//             fout->openGroup(hdf5path).createAttribute(
//                 "t", PredType::IEEE_F64LE, DataSpace(H5S_SCALAR)));
//         attr_t.write(PredType::IEEE_F64LE, &t);
// 
//         const Position3 edge_lengths = space_.edge_lengths();
//         const hsize_t dims[] = {3};
//         const ArrayType lengths_type(PredType::NATIVE_DOUBLE, 1, dims);
//         Attribute attr_lengths(
//             fout->openGroup(hdf5path).createAttribute(
//                 "edge_lengths", lengths_type, DataSpace(H5S_SCALAR)));
//         double lengths[] = {edge_lengths[0], edge_lengths[1], edge_lengths[2]};
//         attr_lengths.write(lengths_type, lengths);
//     }
// 
//     // void save(const std::string& filename)
//     // {
//     //     boost::scoped_ptr<H5::H5File>
//     //         fout(new H5::H5File(filename, H5F_ACC_TRUNC));
// 
//     //     std::ostringstream ost_hdf5path;
//     //     ost_hdf5path << "/" << space_.t();
// 
//     //     boost::scoped_ptr<H5::Group> parent_group(
//     //         new H5::Group(fout->createGroup(ost_hdf5path.str())));
//     //     ost_hdf5path << "/LatticeSpace";
//     //     boost::scoped_ptr<H5::Group>
//     //         group(new H5::Group(parent_group->createGroup(ost_hdf5path.str())));
// 
//     //     save(fout.get(), ost_hdf5path.str());
//     // }
// 
// protected:
// 
//     const space_type& space_;
// };

} // ecell4

#endif /*  __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP */
