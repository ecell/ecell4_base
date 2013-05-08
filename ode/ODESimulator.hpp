#ifndef __ECELL4_ODE_ODE_SIMULATOR_HPP
#define __ECELL4_ODE_ODE_SIMULATOR_HPP

#include <cstring>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/get_mapper_mf.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Simulator.hpp>

#include "ODEWorld.hpp"

#include <H5Cpp.h>
#include <hdf5.h>

#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

namespace ecell4
{

namespace ode
{

class ODESystem
{
public:

    typedef std::vector<double> state_type;

protected:

    typedef utils::get_mapper_mf<
        Species, state_type::size_type>::type species_map_type;

public:

    ODESystem(boost::shared_ptr<NetworkModel> model, const Real& volume)
        : model_(model), volume_(volume)
    {
        initialize();
    }

    void initialize()
    {
        const NetworkModel::species_container_type& species(model_->species());
        state_type::size_type i(0);
        for (NetworkModel::species_container_type::const_iterator
                 it(species.begin()); it != species.end(); ++it)
        {
            index_map_[*it] = i;
            ++i;
        }
    }

    void operator()(const state_type& x, state_type& dxdt, const double& t)
    {
        for (state_type::iterator i(dxdt.begin()); i != dxdt.end(); ++i)
        {
            *i = 0.0;
        }

        const NetworkModel::reaction_rule_container_type&
            reaction_rules(model_->reaction_rules());
        for (NetworkModel::reaction_rule_container_type::const_iterator
                 i(reaction_rules.begin()); i != reaction_rules.end(); ++i)
        {
            double flux((*i).k() * volume_);

            const ReactionRule::reactant_container_type&
                reactants((*i).reactants());
            const ReactionRule::product_container_type&
                products((*i).products());
            for (ReactionRule::reactant_container_type::iterator
                     j(reactants.begin()); j != reactants.end(); ++j)
            {
                flux *= x[index_map_[*j]] / volume_;
            }

            for (ReactionRule::reactant_container_type::iterator
                     j(reactants.begin()); j != reactants.end(); ++j)
            {
                dxdt[index_map_[*j]] -= flux;
            }

            for (ReactionRule::product_container_type::iterator
                     j(products.begin()); j != products.end(); ++j)
            {
                dxdt[index_map_[*j]] += flux;
            }
        }
    }

protected:

    boost::shared_ptr<NetworkModel> model_;
    Real volume_;

    species_map_type index_map_;
};

struct StateAndTimeBackInserter
{
    typedef std::vector<ODESystem::state_type> state_container_type;
    typedef std::vector<double> time_container_type;

    state_container_type& m_states;
    time_container_type& m_times;

    StateAndTimeBackInserter(
        state_container_type& states, time_container_type& times)
        : m_states(states), m_times(times)
    {
        ;
    }

    void operator()(const ODESystem::state_type&x, double t)
    {
        m_states.push_back(x);
        m_times.push_back(t);
    }
};

class ODESimulator
    : public Simulator
{
public:

    ODESimulator(
        boost::shared_ptr<NetworkModel> model,
        boost::shared_ptr<ODEWorld> world)
        : model_(model), world_(world), dt_(0.0), num_steps_(0), is_dirty_(true)
    {
        ;
    }

    void initialize()
    {
        if (!is_dirty_)
        {
            return;
        }

        const NetworkModel::species_container_type& species((*model_).species());
        for (NetworkModel::species_container_type::const_iterator
                 i(species.begin()); i != species.end(); ++i)
        {
            if (!(*world_).has_species(*i))
            {
                (*world_).add_species(*i);
            }
        }

        is_dirty_ = false;
    }

    // SimulatorTraits

    Real t(void) const
    {
        return (*world_).t();
    }

    Integer num_steps(void) const
    {
        return num_steps_;
    }

    Real dt() const
    {
        return dt_;
    }

    void step(void)
    {
        step(next_time());
    }

    bool step(const Real& upto);

    // Optional members

    void set_t(const Real& t)
    {
        (*world_).set_t(t);
    }

    void set_dt(const Real& dt)
    {
        if (dt <= 0)
        {
            throw std::invalid_argument("The step size must be positive.");
        }
        dt_ = dt;
    }

    void save_hdf5_init(std::string filename)
    {
        using namespace H5;
        this->file_ = new H5File(filename, H5F_ACC_TRUNC);
        boost::scoped_ptr<Group>
            group(new Group(this->file_->createGroup("/ODEWorld")));
    }
    void save_hdf5(void)
    {
        using namespace H5;
        if (this->file_ == NULL)
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
            PredType::IEEE_F64LE);

        // Construct Data Set.
        const NetworkModel::species_container_type &species_list =
            this->model_->species();
        boost::scoped_array<species_id_table_struct> species_id_table(
            new species_id_table_struct[species_list.size()]);
        boost::scoped_array<species_num_struct> species_num_table(
            new species_num_struct[species_list.size()]);

        for (unsigned int i(0); i < species_list.size(); i++)
        {
            species_id_table[i].id = i + 1;
            std::strcpy(species_id_table[i].name,
                        species_list[i].name().c_str());

            species_num_table[i].id = i + 1;
            species_num_table[i].num_of_molecules =
                this->world_->num_molecules(species_list[i]);
        }

        const int RANK = 1;
        hsize_t dim[1];
        dim[0] = species_list.size();

        // Create Path.
        std::ostringstream ost_hdf5path;
        boost::scoped_ptr<Group> parent_group(
            new Group(this->file_->openGroup("/ODEWorld")));
        ost_hdf5path << "//ODEWorld/" << this->t();
        boost::scoped_ptr<Group> group(
            new Group(parent_group->createGroup(ost_hdf5path.str())));

        DataSpace space(RANK, dim);
        std::string species_table_path = ost_hdf5path.str() + "/species";
        std::string species_num_path = ost_hdf5path.str() + "/num";
        boost::scoped_ptr<DataSet> dataset_id_table(
            new DataSet(
                this->file_->createDataSet(
                    species_table_path, mtype_id_table_struct, space)));
        boost::scoped_ptr<DataSet> dataset_num_table(
            new DataSet(
                this->file_->createDataSet(
                    species_num_path , mtype_num_struct, space)));

        /* attribute */
        const double t_value = this->t();
        FloatType doubleType(PredType::IEEE_F64LE);

        Attribute attr_id_table = dataset_id_table->createAttribute(
            "t", doubleType, DataSpace(H5S_SCALAR));
        attr_id_table.write(doubleType, &t_value);

        Attribute attr_num_table = dataset_num_table->createAttribute(
            "t", doubleType, DataSpace(H5S_SCALAR));
        attr_num_table.write(doubleType, &t_value);

        /* write */
        dataset_id_table->write(species_id_table.get(), mtype_id_table_struct);
        dataset_num_table->write(species_num_table.get(), mtype_num_struct);
    }

protected:

    typedef struct species_id_table_struct {
            uint32_t id;
            char name[32];
    } species_id_table_struct;

    typedef struct species_num_struct {
            uint32_t id;
            double num_of_molecules;
    } species_num_struct;

    boost::shared_ptr<NetworkModel> model_;
    boost::shared_ptr<ODEWorld> world_;
    Real dt_;
    Integer num_steps_;
    bool is_dirty_;

    H5::H5File *file_;
};

} // ode

} // ecell4

#endif /* __ECELL4_ODE_ODE_SIMULATOR_HPP */
