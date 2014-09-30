#ifndef __ECELL4_SUBVOLUME_SPACE_HDF5_WRITER_HPP
#define __ECELL4_SUBVOLUME_SPACE_HDF5_WRITER_HPP

#include <cstring>
#include <iostream>
#include <sstream>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/multi_array.hpp>
#include <boost/lexical_cast.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"
#include "Species.hpp"
#include "Voxel.hpp"


namespace ecell4
{

struct SubvolumeSpaceHDF5Traits
{
    typedef struct h5_species_struct {
        uint32_t id;
        char serial[32]; // species' serial may exceed the limit
        // double D;
    } h5_species_struct;

    static H5::CompType get_species_comp_type()
    {
        H5::CompType h5_species_comp_type(sizeof(h5_species_struct));
        h5_species_comp_type.insertMember(
            std::string("id"), HOFFSET(h5_species_struct, id),
            H5::PredType::STD_I32LE);
        h5_species_comp_type.insertMember(
            std::string("serial"), HOFFSET(h5_species_struct, serial),
            H5::StrType(H5::PredType::C_S1, 32));
        // h5_species_comp_type.insertMember(
        //     std::string("D"), HOFFSET(h5_species_struct, D),
        //     H5::PredType::STD_I64LE);
        return h5_species_comp_type;
    }
};

template<typename Tspace_>
void save_subvolume_space(const Tspace_& space, H5::Group* root)
{
    typedef SubvolumeSpaceHDF5Traits traits_type;
    typedef typename traits_type::h5_species_struct h5_species_struct;
    // typedef typename traits_type::h5_voxel_struct h5_voxel_struct;

    // typedef std::vector<std::pair<ParticleID, Voxel> >
    //     voxel_container_type;

    const unsigned int num_subvolumes(space.num_subvolumes());
    const std::vector<Species> species(space.list_species());
    boost::multi_array<int64_t, 2>
        h5_num_table(boost::extents[species.size()][num_subvolumes]);
    boost::scoped_array<h5_species_struct>
        h5_species_table(new h5_species_struct[species.size()]);

    for (unsigned int i(0); i < species.size(); ++i)
    {
        const unsigned int sid(i + 1);
        h5_species_table[i].id = sid;
        std::strcpy(h5_species_table[i].serial, species[i].serial().c_str());
        // h5_species_table[i].D = 0.0;
        for (unsigned int j(0); j < num_subvolumes; ++j)
        {
            h5_num_table[i][j] = space.num_molecules_exact(species[i], j);
        }
    }

    const int RANK1 = 2;
    const int RANK2 = 1;

    hsize_t dim1[] = {species.size(), num_subvolumes};
    H5::DataSpace dataspace1(RANK1, dim1);
    boost::scoped_ptr<H5::DataSet> dataset1(new H5::DataSet(
        root->createDataSet(
            "num_molecules", H5::PredType::STD_I64LE, dataspace1)));

    hsize_t dim2[] = {species.size()};
    H5::DataSpace dataspace2(RANK2, dim2);
    boost::scoped_ptr<H5::DataSet> dataset2(new H5::DataSet(
        root->createDataSet(
            "species", traits_type::get_species_comp_type(), dataspace2)));

    dataset1->write(h5_num_table.data(), dataset1->getDataType());
    dataset2->write(h5_species_table.get(), dataset2->getDataType());

    const uint32_t space_type = static_cast<uint32_t>(Space::SUBVOLUME);
    H5::Attribute attr_space_type(
        root->createAttribute(
            "type", H5::PredType::STD_I32LE, H5::DataSpace(H5S_SCALAR)));
    attr_space_type.write(H5::PredType::STD_I32LE, &space_type);

    const double t = space.t();
    H5::Attribute attr_t(
        root->createAttribute(
            "t", H5::PredType::IEEE_F64LE, H5::DataSpace(H5S_SCALAR)));
    attr_t.write(H5::PredType::IEEE_F64LE, &t);

    const Position3 edge_lengths = space.edge_lengths();
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    H5::Attribute attr_lengths(
        root->createAttribute(
            "edge_lengths", lengths_type, H5::DataSpace(H5S_SCALAR)));
    double lengths[] = {edge_lengths[0], edge_lengths[1], edge_lengths[2]};
    attr_lengths.write(lengths_type, lengths);

    const Global matrix_sizes = space.matrix_sizes();
    const H5::ArrayType sizes_type(H5::PredType::STD_I64LE, 1, dims);
    H5::Attribute attr_sizes(
        root->createAttribute(
            "matrix_sizes", sizes_type, H5::DataSpace(H5S_SCALAR)));
    int64_t sizes[] = {matrix_sizes.col, matrix_sizes.row, matrix_sizes.layer};
    attr_sizes.write(sizes_type, sizes);
}

template<typename Tspace_>
void load_subvolume_space(const H5::Group& root, Tspace_* space)
{
    typedef SubvolumeSpaceHDF5Traits traits_type;
    typedef typename traits_type::h5_species_struct h5_species_struct;
    // typedef typename traits_type::h5_voxel_struct h5_voxel_struct;

    double t;
    root.openAttribute("t").read(H5::PredType::IEEE_F64LE, &t);
    space->set_t(t);

    Position3 edge_lengths;
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    root.openAttribute("edge_lengths").read(lengths_type, &edge_lengths);

    int64_t sizes[3];
    const H5::ArrayType sizes_type(H5::PredType::STD_I64LE, 1, dims);
    root.openAttribute("matrix_sizes").read(sizes_type, sizes);
    const Global matrix_sizes(sizes[0], sizes[1], sizes[2]);

    space->cleanup(edge_lengths, matrix_sizes);

    {
        H5::DataSet species_dset(root.openDataSet("species"));
        const unsigned int num_species(
            species_dset.getSpace().getSimpleExtentNpoints());
        boost::scoped_array<h5_species_struct> h5_species_table(
            new h5_species_struct[num_species]);
        species_dset.read(
            h5_species_table.get(), traits_type::get_species_comp_type());
        species_dset.close();

        H5::DataSet num_dset(root.openDataSet("num_molecules"));
        hsize_t dims[2];
        num_dset.getSpace().getSimpleExtentDims(dims);
        assert(num_species == dims[0]);
        const unsigned int num_subvolumes(dims[1]);
        boost::multi_array<int64_t, 2>
            h5_num_table(boost::extents[num_species][num_subvolumes]);
        num_dset.read(
            h5_num_table.data(), H5::PredType::STD_I64LE);
        num_dset.close();

        typedef utils::get_mapper_mf<unsigned int, Species::serial_type>::type
            species_id_map_type;
        species_id_map_type species_id_map;
        for (unsigned int i(0); i < num_species; ++i)
        {
            species_id_map[h5_species_table[i].id] = h5_species_table[i].serial;
        }

        for (unsigned int i(0); i < num_species; ++i)
        {
            const uint32_t sid(i + 1);
            const Species sp(species_id_map[sid]);
            for (unsigned int j(0); j < num_subvolumes; ++j)
            {
                space->add_molecules(sp, h5_num_table[i][j], j);
            }
        }
    }
}

} // ecell4

#endif /*  __ECELL4_SUBVOLUME_SPACE_HDF5_WRITER_HPP */
