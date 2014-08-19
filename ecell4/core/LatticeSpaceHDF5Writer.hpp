#ifndef __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP
#define __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP

#include <cstring>
#include <iostream>
#include <sstream>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/lexical_cast.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"
#include "Species.hpp"
#include "Voxel.hpp"


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
        double radius;
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
            std::string("radius"), HOFFSET(h5_voxel_struct, radius),
            H5::PredType::NATIVE_DOUBLE);
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

        const voxel_container_type& voxels(space.list_voxels_exact(species[i]));
        for (unsigned int j(0); j < voxels.size(); ++j)
        {
            h5_voxel_table[vidx].lot = voxels[j].first.lot();
            h5_voxel_table[vidx].serial = voxels[j].first.serial();
            h5_voxel_table[vidx].sid = sid;
            h5_voxel_table[vidx].coord = voxels[j].second.coordinate();
            h5_voxel_table[vidx].radius = voxels[j].second.radius();
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

    const uint32_t is_periodic = (space.is_periodic()? 1 : 0);
    H5::Attribute attr_is_periodic(
        root->createAttribute(
            "is_periodic", H5::PredType::STD_I32LE, H5::DataSpace(H5S_SCALAR)));
    attr_is_periodic.write(H5::PredType::STD_I32LE, &is_periodic);

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

template<typename Tspace_>
void load_lattice_space(const H5::Group& root, Tspace_* space)
{
    typedef LatticeSpaceHDF5Traits traits_type;
    typedef typename traits_type::h5_species_struct h5_species_struct;
    typedef typename traits_type::h5_voxel_struct h5_voxel_struct;

    double t;
    root.openAttribute("t").read(H5::PredType::IEEE_F64LE, &t);
    space->set_t(t);

    uint32_t is_periodic;
    root.openAttribute("is_periodic").read(
        H5::PredType::STD_I32LE, &is_periodic);

    Position3 edge_lengths;
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    root.openAttribute("edge_lengths").read(lengths_type, &edge_lengths);

    double voxel_radius;
    root.openAttribute("voxel_radius").read(H5::PredType::IEEE_F64LE, &voxel_radius);

    space->cleanup(edge_lengths, voxel_radius, (is_periodic != 0));

    {
        H5::DataSet species_dset(root.openDataSet("species"));
        const unsigned int num_species(
            species_dset.getSpace().getSimpleExtentNpoints());
        boost::scoped_array<h5_species_struct> h5_species_table(
            new h5_species_struct[num_species]);
        species_dset.read(
            h5_species_table.get(), traits_type::get_species_comp_type());
        species_dset.close();

        H5::DataSet voxel_dset(root.openDataSet("voxels"));
        const unsigned int num_voxels(
            voxel_dset.getSpace().getSimpleExtentNpoints());
        boost::scoped_array<h5_voxel_struct> h5_voxel_table(
            new h5_voxel_struct[num_voxels]);
        voxel_dset.read(
            h5_voxel_table.get(), traits_type::get_voxel_comp_type());
        voxel_dset.close();

        typedef utils::get_mapper_mf<unsigned int, Species::serial_type>::type
            species_id_map_type;
        species_id_map_type species_id_map;
        for (unsigned int i(0); i < num_species; ++i)
        {
            species_id_map[h5_species_table[i].id] = h5_species_table[i].serial;
        }

        for (unsigned int i(0); i < num_voxels; ++i)
        {
            space->update_voxel(
                ParticleID(std::make_pair(h5_voxel_table[i].lot, h5_voxel_table[i].serial)),
                Voxel(Species(species_id_map[h5_voxel_table[i].sid]),
                    h5_voxel_table[i].coord, h5_voxel_table[i].radius,
                    h5_voxel_table[i].D));
        }
    }
}

} // ecell4

#endif /*  __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP */
