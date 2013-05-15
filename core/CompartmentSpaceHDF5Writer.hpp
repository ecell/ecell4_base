#ifndef __ECELL4_COMPARTMENT_SPACE_HDF5_WRITER_HPP
#define __ECELL4_COMPARTMENT_SPACE_HDF5_WRITER_HPP

#include <cstring>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "types.hpp"
#include "Species.hpp"


namespace ecell4
{

struct H5DataTypeTraits_uint32_t
{
    typedef uint32_t type;

    const H5::DataType& operator()() const
    {
        return H5::PredType::STD_I32LE;
    }
};

struct H5DataTypeTraits_double
{
    typedef double type;

    const H5::DataType& operator()() const
    {
        return H5::PredType::IEEE_F64LE;
    }
};

template<typename Tspace_, typename Ttraits_ = H5DataTypeTraits_uint32_t>
class CompartmentSpaceHDF5Writer
{
public:

    typedef Tspace_ space_type;
    typedef Ttraits_ traits_type;
    typedef typename Ttraits_::type num_molecules_type;

protected:

    typedef struct species_id_table_struct {
        uint32_t sid;
        char serial[32]; // species' serial may exceed the limit
    } species_id_table_struct;

    typedef struct species_num_struct {
        uint32_t sid;
        num_molecules_type num_molecules;
    } species_num_struct;

public:

    CompartmentSpaceHDF5Writer(const space_type& space)
        : space_(space), traits_()
    {
        ;
    }

    virtual ~CompartmentSpaceHDF5Writer()
    {
        ;
    }

    void save(H5::H5File* fout, const std::string& hdf5path) const
    {
        using namespace H5;

        const std::vector<Species> species_list(space_.list_species());
        const std::vector<Species>::size_type num_species(species_list.size());

        CompType mtype_id_table_struct(sizeof(species_id_table_struct));
        mtype_id_table_struct.insertMember(
            std::string("sid"), HOFFSET(species_id_table_struct, sid),
            PredType::STD_I32LE);
        mtype_id_table_struct.insertMember(
            std::string("serial"), HOFFSET(species_id_table_struct, serial),
            StrType(PredType::C_S1, 32));

        CompType mtype_num_struct(sizeof(species_num_struct));
        mtype_num_struct.insertMember(
            std::string("sid"), HOFFSET(species_num_struct, sid),
            PredType::STD_I32LE);
        mtype_num_struct.insertMember(
            std::string("num_molecules"),
            HOFFSET(species_num_struct, num_molecules),
            traits_());

        boost::scoped_array<species_id_table_struct>
            species_id_table(new species_id_table_struct[num_species]);
        boost::scoped_array<species_num_struct>
            species_num_table(new species_num_struct[num_species]);

        for(unsigned int i(0); i < num_species; ++i)
        {
            species_id_table[i].sid = i + 1;
            std::strcpy(species_id_table[i].serial,
                        species_list[i].serial().c_str());

            species_num_table[i].sid = i + 1;
            species_num_table[i].num_molecules =
                space_.num_molecules(species_list[i]);
        }

        const int RANK = 1;
        hsize_t dim[1];
        dim[0] = num_species;

        DataSpace dataspace(RANK, dim);
        const std::string
            species_table_path(hdf5path + "/species"),
            species_num_path(hdf5path + "/num_molecules");
        boost::scoped_ptr<DataSet> dataset_id_table(
            new DataSet(fout->createDataSet(
                            species_table_path, mtype_id_table_struct,
                            dataspace)));
        boost::scoped_ptr<DataSet> dataset_num_table(
            new DataSet(fout->createDataSet(
                            species_num_path, mtype_num_struct, dataspace)));
        dataset_id_table->write(species_id_table.get(), mtype_id_table_struct);
        dataset_num_table->write(species_num_table.get(), mtype_num_struct);

        // attributes
        const double t = space_.t();
        Attribute attr_t(
            fout->openGroup(hdf5path).createAttribute(
                "t", PredType::IEEE_F64LE, DataSpace(H5S_SCALAR)));
        attr_t.write(PredType::IEEE_F64LE, &t);

        const double volume = space_.volume();
        Attribute attr_volume(
            fout->openGroup(hdf5path).createAttribute(
                "volume", PredType::IEEE_F64LE, DataSpace(H5S_SCALAR)));
        attr_volume.write(PredType::IEEE_F64LE, &volume);
    }

    // void save(const std::string& filename)
    // {
    //     boost::scoped_ptr<H5::H5File>
    //         fout(new H5::H5File(filename, H5F_ACC_TRUNC));

    //     std::ostringstream ost_hdf5path;
    //     ost_hdf5path << "/" << space_.t();

    //     boost::scoped_ptr<H5::Group> parent_group(
    //         new H5::Group(fout->createGroup(ost_hdf5path.str())));
    //     ost_hdf5path << "/CompartmentSpace";
    //     boost::scoped_ptr<H5::Group>
    //         group(new H5::Group(parent_group->createGroup(ost_hdf5path.str())));

    //     save(fout.get(), ost_hdf5path.str());
    // }

protected:

    const space_type& space_;
    const traits_type traits_;
};

} // ecell4

#endif /*  __ECELL4_COMPARTMENT_SPACE_HDF5_WRITER_HPP */
