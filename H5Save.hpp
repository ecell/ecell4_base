
#include <vector>
#include <string>
#include <sstream>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#include <ecell4/core/types.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Species.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

using namespace H5;


namespace ecell4 {

//============================================================
//	get_h5_scalar_data_type_le
//============================================================
template<typename T>
struct get_h5_scalar_data_type_le
{
    void operator()() const
    {
    }
};

#define MAP_H5_SCALAR_TYPE_LE(type, h5type) \
    template<> struct get_h5_scalar_data_type_le<type> \
    { H5::DataType operator()() const { return h5type; } }

MAP_H5_SCALAR_TYPE_LE(char, H5::PredType::C_S1);
MAP_H5_SCALAR_TYPE_LE(uint8_t, H5::PredType::STD_U8LE);
MAP_H5_SCALAR_TYPE_LE(uint16_t, H5::PredType::STD_U16LE);
MAP_H5_SCALAR_TYPE_LE(uint32_t, H5::PredType::STD_U32LE);
MAP_H5_SCALAR_TYPE_LE(uint64_t, H5::PredType::STD_U64LE);
MAP_H5_SCALAR_TYPE_LE(int8_t, H5::PredType::STD_I8LE);
MAP_H5_SCALAR_TYPE_LE(int16_t, H5::PredType::STD_I16LE);
MAP_H5_SCALAR_TYPE_LE(int32_t, H5::PredType::STD_I32LE);
MAP_H5_SCALAR_TYPE_LE(int64_t, H5::PredType::STD_I64LE);
MAP_H5_SCALAR_TYPE_LE(float, H5::PredType::IEEE_F32LE);
MAP_H5_SCALAR_TYPE_LE(double, H5::PredType::IEEE_F64LE);

#undef MAP_H5_SCALAR_TYPE_LE

//============================================================
//
//============================================================
typedef struct species_id_table_struct {
	uint32_t id;
	char name[32];
} species_id_table_struct;

struct mtype_species_id_table_struct 
{
	CompType mtype;
	mtype_species_id_table_struct(void)
	{
		CompType mtype_id_table_struct(sizeof(species_id_table_struct));
		mtype_id_table_struct.insertMember
				(std::string("id"), HOFFSET(species_id_table_struct, id), PredType::STD_I32LE);
		mtype_id_table_struct.insertMember
				(std::string("name"), HOFFSET(species_id_table_struct, name), 
				 StrType(PredType::C_S1, 32));
		this->mtype = mtype_id_table_struct;
	}
	CompType operator()(void)
	{
		return this->mtype;
	}
};

//============================================================
//
//============================================================
template <typename T = double>
struct species_num_struct {
	uint32_t id;
	T num_of_molecules;
};

template<typename T>
struct mtype_species_num_struct 
{
	CompType mtype;
	struct get_h5_scalar_data_type_le<T> typeT;
	//Constructor
	mtype_species_num_struct(void)
	{
		CompType mtype_id_table_struct(sizeof(species_num_struct<T>));
		mtype_id_table_struct.insertMember
				(std::string("id"), HOFFSET(species_num_struct<T>, id), PredType::STD_I32LE);
		mtype_id_table_struct.insertMember
				(std::string("number"), HOFFSET(species_num_struct<T>, num_of_molecules), typeT());
		this->mtype = mtype_id_table_struct;
	}
	CompType operator() (void)
	{
		return this->mtype;
	}
};


//============================================================
//
//============================================================
template <typename WorldType, typename MoleculesNumType>
class ecell4_hdf5_manager 
{
	H5File *file_;

	std::string group_root_;
	mtype_species_id_table_struct table_type_;
	boost::shared_ptr<ecell4::NetworkModel> model_;

	// Template Parameter Dependent.
	boost::shared_ptr<WorldType> world_;
	mtype_species_num_struct<MoleculesNumType> num_type_;

public:
	// Constructor
	ecell4_hdf5_manager(
					std::string filename, 
					boost::shared_ptr<NetworkModel> m, 
					boost::shared_ptr<WorldType> w,
					std::string group_name = std::string("World") )
			:model_(m) , world_(w) 
	{
		using namespace H5;
		this->file_ = new H5File(filename, H5F_ACC_TRUNC);

		// set group path.
		this->group_root_ = std::string("/") + group_name;
		boost::scoped_ptr<Group> group (new Group(this->file_->createGroup(group_root_)));

		this->model_ = m;
		this->world_ = w;
	}

	// Destructor
	~ecell4_hdf5_manager(void)
	{
		delete this->file_;
	}

	void save(void)
	{
		// Construct Data Set.
		const NetworkModel::species_container_type &species_list = this->model_->species();
		boost::scoped_array<species_id_table_struct> species_id_table(new species_id_table_struct[species_list.size() ]);
		boost::scoped_array<species_num_struct<MoleculesNumType> > species_num_table(new species_num_struct<MoleculesNumType> [ species_list.size() ]);
		
		for(unsigned int i(0); i < species_list.size(); i++) {
			species_id_table[i].id = i + 1;
			std::strcpy(species_id_table[i].name, species_list[i].name().c_str());
		
			species_num_table[i].id = i + 1;
			species_num_table[i].num_of_molecules = this->world_->num_molecules( species_list[i] );
		}
		const int RANK = 1;
		hsize_t dim[1];
		dim[0] = species_list.size();
		
		// Create Path.
		std::ostringstream ost_hdf5path;
		boost::scoped_ptr<Group> parent_group (new Group(this->file_->openGroup(this->group_root_)));
		ost_hdf5path << this->group_root_ + "/" << this->world_->t();
		boost::scoped_ptr<Group> group (new Group(parent_group->createGroup( ost_hdf5path.str() )));
		
		DataSpace space(RANK, dim);
		std::string species_table_path = ost_hdf5path.str() + "/species";
		std::string species_num_path = ost_hdf5path.str() + "/num";
		boost::scoped_ptr<DataSet> dataset_id_table( new DataSet(this->file_->createDataSet(species_table_path , this->table_type_(), space)) );
		boost::scoped_ptr<DataSet> dataset_num_table( new DataSet(this->file_->createDataSet(species_num_path , this->num_type_(), space)) );
		
		/* set attribute */
		const double t_value = this->world_->t();
		FloatType doubleType(PredType::IEEE_F64LE);
		
		Attribute attr_id_table = dataset_id_table->createAttribute("t", doubleType, DataSpace(H5S_SCALAR));
		attr_id_table.write(doubleType, &t_value);
		
		Attribute attr_num_table = dataset_num_table->createAttribute("t", doubleType, DataSpace(H5S_SCALAR));
		attr_num_table.write(doubleType, &t_value);
		
		/* write */
		dataset_id_table->write(species_id_table.get(), this->table_type_());
		dataset_num_table->write(species_num_table.get(), this->num_type_());
	}
};


}	// ecell4

