#ifndef __ECELL4_GILLESPIE_COMPARTMENT_SPACE_HDF5_WRITER_HPP
#define __ECELL4_GILLESPIE_COMPARTMENT_SPACE_HDF5_WRITER_HPP

#include <cstring>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/NetworkModel.hpp>


namespace ecell4
{

template<typename Tspace_>
class CompartmentSpaceHDF5Writer
{
public:

    typedef Tspace_ space_type;

protected:

    typedef struct species_id_table_struct {
        uint32_t id;
        char name[32];
    } species_id_table_struct;

    typedef struct species_num_struct {
        uint32_t id;
        uint32_t num_of_molecules;
    } species_num_struct;

public:

    CompartmentSpaceHDF5Writer(const space_type& space)
        : space_(space)
    {
        ;
    }

    virtual ~CompartmentSpaceHDF5Writer()
    {
        if (fout_ != NULL)
        {
            delete fout_;
        }
    }

    void initialize(const std::string& filename)
    {
        using namespace H5;
        fout_ = new H5File(filename, H5F_ACC_TRUNC);
        boost::scoped_ptr<Group> group(
            new Group(fout_->createGroup("/CompartmentSpace")));
    }

    void save()
    {
        using namespace H5;

        if (fout_ == NULL)
        {
            return;
        }

        // Define Data Structure
        CompType mtype_id_table_struct(sizeof(species_id_table_struct));
        mtype_id_table_struct.insertMember(
            std::string("id"), HOFFSET(species_id_table_struct, id),
            PredType::STD_I32LE);
        mtype_id_table_struct.insertMember(
            std::string("name"), HOFFSET(species_id_table_struct, name),
            StrType(PredType::C_S1, 32));

        CompType mtype_num_struct(sizeof(species_num_struct));
        mtype_num_struct.insertMember(
            std::string("id"), HOFFSET(species_num_struct, id),
            PredType::STD_I32LE);
        mtype_num_struct.insertMember(
            std::string("number"), HOFFSET(species_num_struct, num_of_molecules),
            PredType::STD_I32LE);

        // Construct Data Set.
        const NetworkModel::species_container_type&
            species_list(space_.list_species());

        boost::scoped_array<species_id_table_struct>
            species_id_table(new species_id_table_struct[species_list.size()]);
        boost::scoped_array<species_num_struct>
            species_num_table(new species_num_struct[species_list.size()]);

        for(unsigned int i(0); i < species_list.size(); ++i)
        {
            species_id_table[i].id = i + 1;
            std::strcpy(species_id_table[i].name, species_list[i].name().c_str());

            species_num_table[i].id = i + 1;
            species_num_table[i].num_of_molecules =
                space_.num_molecules(species_list[i]);
        }

        const int RANK = 1;
        hsize_t dim[1];
        dim[0] = species_list.size();

        // Create Path.
        std::ostringstream ost_hdf5path;
        boost::scoped_ptr<Group>
            parent_group(new Group(fout_->openGroup("/CompartmentSpace")));
        ost_hdf5path << "/CompartmentSpace/" << space_.t();
        boost::scoped_ptr<Group>
            group(new Group(parent_group->createGroup(ost_hdf5path.str())));

        DataSpace space(RANK, dim);
        std::string species_table_path = ost_hdf5path.str() + "/species";
        std::string species_num_path = ost_hdf5path.str() + "/num";
        boost::scoped_ptr<DataSet> dataset_id_table(
            new DataSet(fout_->createDataSet(
                            species_table_path, mtype_id_table_struct, space)));
        boost::scoped_ptr<DataSet> dataset_num_table(
            new DataSet(fout_->createDataSet(
                            species_num_path, mtype_num_struct, space)));

        // attribute
        const double t_value = space_.t();
        FloatType doubleType(PredType::IEEE_F64LE);

        Attribute attr_id_table(
            dataset_id_table->createAttribute(
                "t", doubleType, DataSpace(H5S_SCALAR)));
        attr_id_table.write(doubleType, &t_value);

        Attribute attr_num_table(
            dataset_num_table->createAttribute(
                "t", doubleType, DataSpace(H5S_SCALAR)));
        attr_num_table.write(doubleType, &t_value);

        // write
        dataset_id_table->write(species_id_table.get(), mtype_id_table_struct);
        dataset_num_table->write(species_num_table.get(), mtype_num_struct);
    }

protected:

    const space_type& space_;
    H5::H5File *fout_;
};

} // ecell4

#endif /*  __ECELL4_GILLESPIE_COMPARTMENT_SPACE_HDF5_WRITER_HPP */
