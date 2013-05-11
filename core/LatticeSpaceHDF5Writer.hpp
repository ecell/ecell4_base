#ifndef __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP
#define __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP

#include <cstring>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"
#include "Species.hpp"
#include "Voxel.hpp"


namespace ecell4
{

template<typename Tspace_>
class LatticeSpaceHDF5Writer
{
public:

    typedef Tspace_ space_type;

protected:

    struct h5_species_struct
    {
        uint32_t id;
        char serial[32]; // species' serial may exceed the limit
    };

    struct h5_voxel_struct
    {
        uint32_t sid;
    };

public:

    LatticeSpaceHDF5Writer(const space_type& space)
        : space_(space)
    {
        ;
    }

    virtual ~LatticeSpaceHDF5Writer()
    {
        ;
    }

    void save(H5::H5File* fout, const std::string& hdf5path) const
    {
        using namespace H5;

        const Integer num_voxels(space_.num_voxels());

        std::vector<Species> species;
        typedef utils::get_mapper_mf<Species::serial_type, unsigned int>::type
        species_id_map_type;
        species_id_map_type species_id_map;

        boost::scoped_array<h5_voxel_struct>
            h5_voxel_table(new h5_voxel_struct[num_voxels]);
        for (unsigned int i(0); i < num_voxels; ++i)
        {
            const Voxel voxel(space_.get_voxel(i));

            species_id_map_type::const_iterator
                it(species_id_map.find(voxel.species.serial()));
            if (it == species_id_map.end())
            {
                species.push_back(voxel.species);
                it = species_id_map.insert(
                    std::make_pair(voxel.species.serial(),
                                   species.size())).first;
            }

            h5_voxel_table[i].sid = (*it).second;
        }

        boost::scoped_array<h5_species_struct>
            h5_species_table(new h5_species_struct[species.size()]);
        for (unsigned int i(0); i < species.size(); ++i)
        {
            h5_species_table[i].id = i + 1;
            std::strncpy(h5_species_table[i].serial,
                         species[i].serial().c_str(),
                         sizeof(h5_species_table[i].serial));
            // std::strcpy(h5_species_table[i].serial,
            //             species[i].serial().c_str());
        }

        CompType h5_voxel_comp_type(sizeof(h5_voxel_struct));
        h5_voxel_comp_type.insertMember(
            std::string("sid"), HOFFSET(h5_voxel_struct, sid),
            PredType::STD_I32LE);

        CompType h5_species_comp_type(sizeof(h5_species_struct));
        h5_species_comp_type.insertMember(
            std::string("id"), HOFFSET(h5_species_struct, id),
            PredType::STD_I32LE);
        h5_species_comp_type.insertMember(
            std::string("serial"), HOFFSET(h5_species_struct, serial),
            StrType(PredType::C_S1, 32));

        const int RANK = 1;
        hsize_t dim1[] = {num_voxels};
        DataSpace dataspace1(RANK, dim1);
        boost::scoped_ptr<DataSet> dataset1(
            new DataSet(fout->createDataSet(
                            hdf5path + "/voxels", h5_voxel_comp_type,
                            dataspace1)));

        hsize_t dim2[] = {species.size()};
        DataSpace dataspace2(RANK, dim2);
        boost::scoped_ptr<DataSet> dataset2(
            new DataSet(fout->createDataSet(
                            hdf5path + "/species" , h5_species_comp_type,
                            dataspace2)));
        dataset1->write(h5_voxel_table.get(), h5_voxel_comp_type);
        dataset2->write(h5_species_table.get(), h5_species_comp_type);

        const double t = space_.t();
        Attribute attr_t(
            fout->openGroup(hdf5path).createAttribute(
                "t", PredType::IEEE_F64LE, DataSpace(H5S_SCALAR)));
        attr_t.write(PredType::IEEE_F64LE, &t);

        const double voxel_radius = space_.voxel_radius();
        Attribute attr_radius(
            fout->openGroup(hdf5path).createAttribute(
                "voxel_radius", PredType::IEEE_F64LE, DataSpace(H5S_SCALAR)));
        attr_radius.write(PredType::IEEE_F64LE, &voxel_radius);

        const Position3 edge_lengths = space_.edge_lengths();
        const hsize_t dims[] = {3};
        const ArrayType lengths_type(PredType::NATIVE_DOUBLE, 1, dims);
        Attribute attr_lengths(
            fout->openGroup(hdf5path).createAttribute(
                "edge_lengths", lengths_type, DataSpace(H5S_SCALAR)));
        double lengths[] = {edge_lengths[0], edge_lengths[1], edge_lengths[2]};
        attr_lengths.write(lengths_type, lengths);
    }

    // void save(const std::string& filename)
    // {
    //     boost::scoped_ptr<H5::H5File>
    //         fout(new H5::H5File(filename, H5F_ACC_TRUNC));

    //     std::ostringstream ost_hdf5path;
    //     ost_hdf5path << "/" << space_.t();

    //     boost::scoped_ptr<H5::Group> parent_group(
    //         new H5::Group(fout->createGroup(ost_hdf5path.str())));
    //     ost_hdf5path << "/LatticeSpace";
    //     boost::scoped_ptr<H5::Group>
    //         group(new H5::Group(parent_group->createGroup(ost_hdf5path.str())));

    //     save(fout.get(), ost_hdf5path.str());
    // }

protected:

    const space_type& space_;
};

} // ecell4

#endif /*  __ECELL4_LATTICE_SPACE_HDF5_WRITER_HPP */
