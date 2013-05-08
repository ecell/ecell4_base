#ifndef __ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP
#define __ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Simulator.hpp>

#include <hdf5.h>
#include <H5Cpp.h>

#include "GillespieWorld.hpp"


namespace ecell4
{

namespace gillespie
{

class GillespieSimulator
    : public Simulator
{
public:

    GillespieSimulator(
        boost::shared_ptr<NetworkModel> model,
        boost::shared_ptr<GillespieWorld> world)
        : model_(model), world_(world), num_steps_(0)
    {
        this->initialize();

        // About Hdf5
        this->file_ = NULL;
    }

    ~GillespieSimulator(void)
    {
        if (this->file_ != NULL)
        {
            delete this->file_;
        }
    }

    // SimulatorTraits

    Real t(void) const;
    Real dt(void) const;
    Integer num_steps(void) const;

    void step(void) ;
    bool step(const Real & upto);

    // Optional members

    void set_t(const Real &t);

    /**
     * recalculate reaction propensities and draw the next time.
     */
    void initialize(void);

    // About Hdf5
    void save_hdf5_init(std::string filename);
    void save_hdf5(void);

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return (*world_).rng();
    }

protected:

    void draw_next_reaction(void);

protected:

    boost::shared_ptr<NetworkModel> model_;
    boost::shared_ptr<GillespieWorld> world_;
    Integer num_steps_;

    Real dt_;
    int next_reaction_num_; // the index of the next reaction.

    // About Hdf5
    H5::H5File *file_;

    typedef struct species_id_table_struct {
        uint32_t id;
        char name[32];
    } species_id_table_struct;

    typedef struct species_num_struct {
        uint32_t id;
        uint32_t num_of_molecules;
    } species_num_struct;
};

}

} // ecell4

#endif /* __ECELL4_GILLESPIE_GILLESPIE_SIMULATOR_HPP */
