#ifndef ECELL4_COMPARTMENT_SPACE_HDF5_WRITER_HPP
#define ECELL4_COMPARTMENT_SPACE_HDF5_WRITER_HPP

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

template<typename Tspace_, typename Tdata_>
struct CompartmentSpaceHDF5TraitsBase
{
    typedef Tspace_ space_type;
    typedef Tdata_ num_molecules_traits_type;
    typedef typename num_molecules_traits_type::type num_molecules_type;

    typedef struct species_id_table_struct {
        uint32_t sid;
        char serial[32]; // species' serial may exceed the limit
    } species_id_table_struct;

    static H5::CompType get_species_id_table_struct_memtype()
    {
        H5::CompType mtype_id_table_struct(sizeof(species_id_table_struct));
        // const H5std_string name1("sid");
        // const H5std_string name2("serial");
        // mtype_id_table_struct.insertMember(
        //     name1, HOFFSET(species_id_table_struct, sid),
        //     H5::PredType::STD_I32LE);
        // mtype_id_table_struct.insertMember(
        //     name2, HOFFSET(species_id_table_struct, serial),
        //     H5::StrType(H5::PredType::C_S1, 32));
#define INSERT_MEMBER(member, type) \
        H5Tinsert(mtype_id_table_struct.getId(), #member,\
                HOFFSET(species_id_table_struct, member), type.getId())
        INSERT_MEMBER(sid, H5::PredType::STD_I32LE);
        INSERT_MEMBER(serial, H5::StrType(H5::PredType::C_S1, 32));
#undef INSERT_MEMBER
        return mtype_id_table_struct;
    }

    typedef struct species_num_struct {
        uint32_t sid;
        num_molecules_type num_molecules;
    } species_num_struct;

    static H5::CompType get_species_num_struct_memtype()
    {
        H5::CompType mtype_num_struct(sizeof(species_num_struct));
        // const H5std_string name1("sid");
        // const H5std_string name2("num_molecules");
        // mtype_num_struct.insertMember(
        //     name1, HOFFSET(species_num_struct, sid),
        //     H5::PredType::STD_I32LE);
        // mtype_num_struct.insertMember(
        //     name2,
        //     HOFFSET(species_num_struct, num_molecules),
        //     num_molecules_traits_type::get());
#define INSERT_MEMBER(member, type) \
        H5Tinsert(mtype_num_struct.getId(), #member,\
                HOFFSET(species_num_struct, member), type.getId())
        INSERT_MEMBER(sid, H5::PredType::STD_I32LE);
        INSERT_MEMBER(num_molecules, num_molecules_traits_type::get());
#undef INSERT_MEMBER
        return mtype_num_struct;
    }

    virtual num_molecules_type getter(
        const space_type& space, const Species& sp) const = 0;
    virtual void setter(
        space_type& space, const Species& sp, const num_molecules_type& value) const = 0;
};

template<typename Tspace_>
struct CompartmentSpaceHDF5Traits
    : public CompartmentSpaceHDF5TraitsBase<Tspace_, H5DataTypeTraits_uint32_t>
{
    typedef CompartmentSpaceHDF5TraitsBase<Tspace_, H5DataTypeTraits_uint32_t> base_type;
    typedef typename base_type::num_molecules_type num_molecules_type;
    typedef typename base_type::space_type space_type;

    num_molecules_type getter(const space_type& space, const Species& sp) const
    {
        return space.num_molecules_exact(sp);
    }

    void setter(
        Tspace_& space, const Species& sp, const num_molecules_type& value) const
    {
        space.add_molecules(sp, value);
    }
};

// template<typename Tspace_, typename Tdata_>
template<typename Ttraits_>
void save_compartment_space(const typename Ttraits_::space_type& space, H5::Group* root)
{
    // typedef CompartmentSpaceHDF5Traits<Tdata_> traits_type;
    typedef Ttraits_ traits_type;
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

    const Real3 edge_lengths = space.edge_lengths();
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    H5::Attribute attr_lengths(
        root->createAttribute(
            "edge_lengths", lengths_type, H5::DataSpace(H5S_SCALAR)));
    double lengths[] = {edge_lengths[0], edge_lengths[1], edge_lengths[2]};
    attr_lengths.write(lengths_type, lengths);
}

// template<typename Tspace_, typename Tdata_>
template<typename Ttraits_>
void load_compartment_space(const H5::Group& root, typename Ttraits_::space_type* space)
{
    // typedef CompartmentSpaceHDF5Traits<Tdata_> traits_type;
    typedef Ttraits_ traits_type;
    typedef typename traits_type::num_molecules_type num_molecules_type;
    typedef typename traits_type::species_id_table_struct species_id_table_struct;
    typedef typename traits_type::species_num_struct species_num_struct;

    Real3 edge_lengths;
    const hsize_t dims[] = {3};
    const H5::ArrayType lengths_type(H5::PredType::NATIVE_DOUBLE, 1, dims);
    root.openAttribute("edge_lengths").read(lengths_type, &edge_lengths);
    space->reset(edge_lengths);

    double t;
    root.openAttribute("t").read(H5DataTypeTraits_double::get(), &t);
    space->set_t(t);

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
}

} // ecell4

#endif /*  ECELL4_COMPARTMENT_SPACE_HDF5_WRITER_HPP */
