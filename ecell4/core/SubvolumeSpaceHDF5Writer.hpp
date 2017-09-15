#ifndef ECELL4_SUBVOLUME_SPACE_HDF5_WRITER_HPP
#define ECELL4_SUBVOLUME_SPACE_HDF5_WRITER_HPP

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
        double D;
        char loc[32]; // species' loc may exceed the limit
    } h5_species_struct;

    static H5::CompType get_species_comp_type()
    {
        H5::CompType h5_species_comp_type(sizeof(h5_species_struct));
#define INSERT_MEMBER(member, type) \
        H5Tinsert(h5_species_comp_type.getId(), #member,\
                HOFFSET(h5_species_struct, member), type.getId())
        INSERT_MEMBER(id, H5::PredType::STD_I32LE);
        INSERT_MEMBER(serial, H5::StrType(H5::PredType::C_S1, 32));
        INSERT_MEMBER(D, H5::PredType::NATIVE_DOUBLE);
        INSERT_MEMBER(loc, H5::StrType(H5::PredType::C_S1, 32));
#undef INSERT_MEMBER
        // h5_species_comp_type.insertMember(
        //     std::string("id"), HOFFSET(h5_species_struct, id),
        //     H5::PredType::STD_I32LE);
        // h5_species_comp_type.insertMember(
        //     std::string("serial"), HOFFSET(h5_species_struct, serial),
        //     H5::StrType(H5::PredType::C_S1, 32));
        // h5_species_comp_type.insertMember(
        //     std::string("D"), HOFFSET(h5_species_struct, D),
        //     H5::PredType::STD_I64LE);  //XXX: -> NATIVE_DOUBLE?
        // h5_species_comp_type.insertMember(
        //     std::string("loc"), HOFFSET(h5_species_struct, loc),
        //     H5::StrType(H5::PredType::C_S1, 32));
        return h5_species_comp_type;
    }

    typedef struct h5_structures_struct {
        uint32_t id;
        char serial[32]; // structures' serial may exceed the limit
        // uint32_t dimension;
    } h5_structures_struct;

    static H5::CompType get_structures_comp_type()
    {
        H5::CompType h5_structures_comp_type(sizeof(h5_structures_struct));
#define INSERT_MEMBER(member, type) \
        H5Tinsert(h5_structures_comp_type.getId(), #member,\
                HOFFSET(h5_structures_struct, member), type.getId())
        INSERT_MEMBER(id, H5::PredType::STD_I32LE);
        INSERT_MEMBER(serial, H5::StrType(H5::PredType::C_S1, 32));
        // INSERT_MEMBER(dimension, H5::PredType::STD_I32LE);
#undef INSERT_MEMBER
        // h5_structures_comp_type.insertMember(
        //     std::string("id"), HOFFSET(h5_structures_struct, id),
        //     H5::PredType::STD_I32LE);
        // h5_structures_comp_type.insertMember(
        //     std::string("serial"), HOFFSET(h5_structures_struct, serial),
        //     H5::StrType(H5::PredType::C_S1, 32));
        // // h5_species_comp_type.insertMember(
        // //     std::string("dimension"), HOFFSET(h5_structure_struct, dimension),
        // //     H5::PredType::STD_I32LE);
        return h5_structures_comp_type;
    }
};

template<typename Tspace_>
void save_subvolume_space(const Tspace_& space, H5::Group* root)
{
    typedef SubvolumeSpaceHDF5Traits traits_type;
    typedef typename traits_type::h5_species_struct h5_species_struct;
    // typedef typename traits_type::h5_voxel_struct h5_voxel_struct;
    typedef typename traits_type::h5_structures_struct h5_structures_struct;

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
        const boost::shared_ptr<typename Tspace_::PoolBase>&
            pool = space.get_pool(species[i]);
        h5_species_table[i].D = pool->D();
        std::strcpy(h5_species_table[i].loc, pool->loc().c_str());

        for (unsigned int j(0); j < num_subvolumes; ++j)
        {
            h5_num_table[i][j] = space.num_molecules_exact(species[i], j);
        }
    }

    const std::vector<Species::serial_type> structures(space.list_structures());
    boost::multi_array<double, 2>
        h5_stcoordinate_table(boost::extents[structures.size()][num_subvolumes]);
    boost::scoped_array<h5_structures_struct>
        h5_structures_table(new h5_structures_struct[structures.size()]);
    for (unsigned int i(0); i < structures.size(); ++i)
    {
        const unsigned int sid(i + 1);
        h5_structures_table[i].id = sid;
        std::strcpy(h5_structures_table[i].serial, structures[i].c_str());
        // h5_structures_table[i].dimension = space.get_dimension(structures[i]);
        for (unsigned int j(0); j < num_subvolumes; ++j)
        {
            // const bool exist = space.check_structure(structures[i], j);
            // h5_stcoordinate_table[i][j] = (exist ? 1 : 0);
            h5_stcoordinate_table[i][j] = space.get_occupancy(structures[i], j);
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

    hsize_t dim3[] = {structures.size(), num_subvolumes};
    H5::DataSpace dataspace3(RANK1, dim3);
    boost::scoped_ptr<H5::DataSet> dataset3(new H5::DataSet(
        root->createDataSet(
            "stcoordinates", H5::PredType::IEEE_F64LE, dataspace3)));

    hsize_t dim4[] = {structures.size()};
    H5::DataSpace dataspace4(RANK2, dim4);
    boost::scoped_ptr<H5::DataSet> dataset4(new H5::DataSet(
        root->createDataSet(
            "structures", traits_type::get_structures_comp_type(), dataspace4)));

    dataset1->write(h5_num_table.data(), dataset1->getDataType());
    dataset2->write(h5_species_table.get(), dataset2->getDataType());
    dataset3->write(h5_stcoordinate_table.data(), dataset3->getDataType());
    dataset4->write(h5_structures_table.get(), dataset4->getDataType());

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

    const Real3 edge_lengths = space.edge_lengths();
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    H5::Attribute attr_lengths(
        root->createAttribute(
            "edge_lengths", lengths_type, H5::DataSpace(H5S_SCALAR)));
    double lengths[] = {edge_lengths[0], edge_lengths[1], edge_lengths[2]};
    attr_lengths.write(lengths_type, lengths);

    const Integer3 matrix_sizes = space.matrix_sizes();
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
    typedef typename traits_type::h5_structures_struct h5_structures_struct;
    // typedef typename traits_type::h5_voxel_struct h5_voxel_struct;

    Real3 edge_lengths;
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    root.openAttribute("edge_lengths").read(lengths_type, &edge_lengths);

    int64_t sizes[3];
    const H5::ArrayType sizes_type(H5::PredType::STD_I64LE, 1, dims);
    root.openAttribute("matrix_sizes").read(sizes_type, sizes);
    const Integer3 matrix_sizes(sizes[0], sizes[1], sizes[2]);

    space->reset(edge_lengths, matrix_sizes);

    double t;
    root.openAttribute("t").read(H5::PredType::IEEE_F64LE, &t);
    space->set_t(t);

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

        typedef utils::get_mapper_mf<unsigned int, unsigned int>::type
            species_id_map_type;
        species_id_map_type species_id_map;
        for (unsigned int i(0); i < num_species; ++i)
        {
            // species_id_map[h5_species_table[i].id] = h5_species_table[i].serial;
            species_id_map[h5_species_table[i].id] = i;
        }

        for (unsigned int i(0); i < num_species; ++i)
        {
            const uint32_t sid(i + 1);
            const unsigned int k(species_id_map[sid]);

            const Species sp(h5_species_table[k].serial);
            const Real D(h5_species_table[k].D);
            const Species::serial_type loc(h5_species_table[k].loc);
            space->reserve_pool(sp, D, loc);

            for (unsigned int j(0); j < num_subvolumes; ++j)
            {
                space->add_molecules(sp, h5_num_table[i][j], j);
            }
        }
    }

    {
        H5::DataSet structures_dset(root.openDataSet("structures"));
        const unsigned int num_structures(
            structures_dset.getSpace().getSimpleExtentNpoints());
        boost::scoped_array<h5_structures_struct> h5_structures_table(
            new h5_structures_struct[num_structures]);
        structures_dset.read(
            h5_structures_table.get(), traits_type::get_structures_comp_type());
        structures_dset.close();

        H5::DataSet stcoordinate_dset(root.openDataSet("stcoordinates"));
        hsize_t dims[2];
        stcoordinate_dset.getSpace().getSimpleExtentDims(dims);
        assert(num_structures == dims[0]);
        const unsigned int num_subvolumes(dims[1]);
        boost::multi_array<double, 2>
            h5_stcoordinate_table(boost::extents[num_structures][num_subvolumes]);
        stcoordinate_dset.read(
            h5_stcoordinate_table.data(), H5::PredType::IEEE_F64LE);
        stcoordinate_dset.close();

        typedef utils::get_mapper_mf<unsigned int, Species::serial_type>::type
            structures_id_map_type;
        structures_id_map_type structures_id_map;
        for (unsigned int i(0); i < num_structures; ++i)
        {
            structures_id_map[h5_structures_table[i].id] = h5_structures_table[i].serial;
        }

        for (unsigned int i(0); i < num_structures; ++i)
        {
            const uint32_t sid(i + 1);
            const Species::serial_type serial(structures_id_map[sid]);
            for (unsigned int j(0); j < num_subvolumes; ++j)
            {
                space->update_structure(serial, j, h5_stcoordinate_table[i][j]);
            }
        }
    }
}

} // ecell4

#endif /*  ECELL4_SUBVOLUME_SPACE_HDF5_WRITER_HPP */
