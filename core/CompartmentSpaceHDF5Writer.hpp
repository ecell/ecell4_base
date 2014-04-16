#ifndef __ECELL4_COMPARTMENT_SPACE_HDF5_WRITER_HPP
#define __ECELL4_COMPARTMENT_SPACE_HDF5_WRITER_HPP

#include <cstring>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "Species.hpp"


namespace ecell4
{

struct H5DataTypeTraits_uint32_t
{
    typedef uint32_t type;

    static const H5::DataType& get()
    {
        return H5::PredType::STD_I32LE;
    }
};

struct H5DataTypeTraits_double
{
    typedef double type;

    static const H5::DataType& get()
    {
        return H5::PredType::IEEE_F64LE;
    }
};

template<typename Tdata_>
struct CompartmentSpaceHDF5Traits
{
    typedef Tdata_ num_molecules_traits_type;
    typedef typename num_molecules_traits_type::type num_molecules_type;

    typedef struct species_id_table_struct {
        uint32_t sid;
        char serial[32]; // species' serial may exceed the limit
    } species_id_table_struct;

    static H5::CompType get_species_id_table_struct_memtype()
    {
        H5::CompType mtype_id_table_struct(sizeof(species_id_table_struct));
        mtype_id_table_struct.insertMember(
            std::string("sid"), HOFFSET(species_id_table_struct, sid),
            H5::PredType::STD_I32LE);
        mtype_id_table_struct.insertMember(
            std::string("serial"), HOFFSET(species_id_table_struct, serial),
            H5::StrType(H5::PredType::C_S1, 32));
        return mtype_id_table_struct;
    }

    typedef struct species_num_struct {
        uint32_t sid;
        num_molecules_type num_molecules;
    } species_num_struct;

    static H5::CompType get_species_num_struct_memtype()
    {
        H5::CompType mtype_num_struct(sizeof(species_num_struct));
        mtype_num_struct.insertMember(
            std::string("sid"), HOFFSET(species_num_struct, sid),
            H5::PredType::STD_I32LE);
        mtype_num_struct.insertMember(
            std::string("num_molecules"),
            HOFFSET(species_num_struct, num_molecules),
            num_molecules_traits_type::get());
        return mtype_num_struct;
    }
};

// template<typename Tspace_, typename Tdata_ = H5DataTypeTraits_uint32_t>
template<typename Tspace_, typename Tdata_>
void save_compartment_space(const Tspace_& space, H5::Group* root)
{
    typedef CompartmentSpaceHDF5Traits<Tdata_> traits_type;
    // typedef typename traits_type::num_molecules_type num_molecules_type;
    typedef typename traits_type::species_id_table_struct species_id_table_struct;
    typedef typename traits_type::species_num_struct species_num_struct;

    // attributes
    const uint32_t space_type = static_cast<uint32_t>(Space::COMPARTMENT);
    H5::Attribute attr_space_type(
        root->createAttribute(
            "type", H5::PredType::STD_I32LE, H5::DataSpace(H5S_SCALAR)));
    attr_space_type.write(H5::PredType::STD_I32LE, &space_type);

    const double t(space.t());
    H5::Attribute attr_t(root->createAttribute(
        "t", H5DataTypeTraits_double::get(), H5::DataSpace(H5S_SCALAR)));
    attr_t.write(attr_t.getDataType(), &t);

    const double volume(space.volume());
    H5::Attribute attr_volume(root->createAttribute(
        "volume", H5DataTypeTraits_double::get(), H5::DataSpace(H5S_SCALAR)));
    attr_volume.write(attr_volume.getDataType(), &volume);

    const std::vector<Species> species_list(space.list_species());
    const std::vector<Species>::size_type num_species(species_list.size());

    boost::scoped_array<species_id_table_struct>
        species_id_table(new species_id_table_struct[num_species]);
    boost::scoped_array<species_num_struct>
        species_num_table(new species_num_struct[num_species]);

    for(unsigned int i(0); i < num_species; ++i)
    {
        species_id_table[i].sid = i + 1;
        std::strcpy(
            species_id_table[i].serial, species_list[i].serial().c_str());

        species_num_table[i].sid = i + 1;
        species_num_table[i].num_molecules =
            space.num_molecules(species_list[i]);
    }

    const int RANK = 1;
    hsize_t dim[1];
    dim[0] = num_species;
    H5::DataSpace dataspace(RANK, dim);

    boost::scoped_ptr<H5::DataSet> dataset_id_table(new H5::DataSet(
        root->createDataSet(
            "species", traits_type::get_species_id_table_struct_memtype(),
            dataspace)));
    boost::scoped_ptr<H5::DataSet> dataset_num_table(new H5::DataSet(
        root->createDataSet(
            "num_molecules", traits_type::get_species_num_struct_memtype(),
            dataspace)));
    dataset_id_table->write(
        species_id_table.get(), dataset_id_table->getDataType());
    dataset_num_table->write(
        species_num_table.get(), dataset_num_table->getDataType());
}

// template<typename Tspace_, typename Tdata_ = H5DataTypeTraits_uint32_t>
template<typename Tspace_, typename Tdata_>
void load_compartment_space(const H5::Group& root, Tspace_* space)
{
    typedef CompartmentSpaceHDF5Traits<Tdata_> traits_type;
    typedef typename traits_type::num_molecules_type num_molecules_type;
    typedef typename traits_type::species_id_table_struct species_id_table_struct;
    typedef typename traits_type::species_num_struct species_num_struct;

    {
        H5::DataSet species_dset(root.openDataSet("species"));
        const unsigned int num_species(
            species_dset.getSpace().getSimpleExtentNpoints());
        boost::scoped_array<species_id_table_struct> species_id_table(
            new species_id_table_struct[num_species]);
        species_dset.read(
            species_id_table.get(),
            traits_type::get_species_id_table_struct_memtype());
        species_dset.close();

        H5::DataSet num_dset(root.openDataSet("num_molecules"));
        boost::scoped_array<species_num_struct> species_num_table(
            new species_num_struct[num_species]);
        num_dset.read(
            species_num_table.get(),
            traits_type::get_species_num_struct_memtype());
        num_dset.close();

        typename utils::get_mapper_mf<uint32_t, num_molecules_type>::type
            num_molecules_cache;
        for (unsigned int i(0); i < num_species; ++i)
        {
            num_molecules_cache[species_num_table[i].sid]
                = species_num_table[i].num_molecules;
        }
        for (unsigned int i(0); i < num_species; ++i)
        {
            space->add_molecules(
                Species(species_id_table[i].serial),
                num_molecules_cache[species_id_table[i].sid]);
        }
    }

    double t;
    root.openAttribute("t").read(H5DataTypeTraits_double::get(), &t);
    space->set_t(t);

    double volume;
    root.openAttribute("volume").read(H5DataTypeTraits_double::get(), &volume);
    space->set_volume(volume);
}

} // ecell4

#endif /*  __ECELL4_COMPARTMENT_SPACE_HDF5_WRITER_HPP */
