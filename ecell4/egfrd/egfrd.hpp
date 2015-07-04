#ifndef __ECELL4_EGFRD_EGFRD_HPP
#define __ECELL4_EGFRD_EGFRD_HPP

#include <ecell4/core/types.hpp>
#include <ecell4/core/SimulatorFactory.hpp>
#include "World.hpp"
#include "EGFRDSimulator.hpp"
#include "BDSimulator.hpp"

namespace ecell4
{

namespace egfrd
{

typedef ::World< ::CyclicWorldTraits<Real> > EGFRDWorld;
typedef EGFRDWorld::molecule_info_type MoleculeInfo;
typedef ::EGFRDSimulator< ::EGFRDSimulatorTraitsBase<EGFRDWorld> > EGFRDSimulator;
typedef ::BDSimulator< ::BDSimulatorTraitsBase<EGFRDWorld> > BDSimulator;

class EGFRDFactory
    : public SimulatorFactory<EGFRDWorld, EGFRDSimulator>
{
public:

    typedef SimulatorFactory<EGFRDWorld, EGFRDSimulator> base_type;

protected:

    typedef EGFRDWorld::matrix_sizes_type matrix_sizes_type;

public:

    EGFRDFactory(
        Integer dissociation_retry_moves=-1, Real bd_dt_factor=-1,
        Real user_max_shell_size=-1)
        : base_type(), matrix_sizes_(0, 0, 0), rng_(),
        num_retries_(dissociation_retry_moves),
        bd_dt_factor_(bd_dt_factor),
        user_max_shell_size_(user_max_shell_size)
    {
        ; // do nothing
    }

    EGFRDFactory(const matrix_sizes_type& matrix_sizes,
        Integer dissociation_retry_moves=-1, Real bd_dt_factor=-1,
        Real user_max_shell_size=-1)
        : base_type(), matrix_sizes_(matrix_sizes), rng_(),
        num_retries_(dissociation_retry_moves),
        bd_dt_factor_(bd_dt_factor),
        user_max_shell_size_(user_max_shell_size)
    {
        ; // do nothing
    }

    EGFRDFactory(const matrix_sizes_type& matrix_sizes,
        const boost::shared_ptr<RandomNumberGenerator>& rng,
        Integer dissociation_retry_moves=-1, Real bd_dt_factor=-1,
        Real user_max_shell_size=-1)
        : base_type(), matrix_sizes_(matrix_sizes), rng_(rng),
        num_retries_(dissociation_retry_moves),
        bd_dt_factor_(bd_dt_factor),
        user_max_shell_size_(user_max_shell_size)
    {
        ; // do nothing
    }

    virtual ~EGFRDFactory()
    {
        ; // do nothing
    }

    virtual EGFRDWorld* create_world(const std::string filename) const
    {
        return new EGFRDWorld(filename);
    }

    virtual EGFRDWorld* create_world(
        const Real3& edge_lengths = Real3(1, 1, 1)) const
    {
        if (rng_)
        {
            return new EGFRDWorld(edge_lengths, matrix_sizes_, rng_);
        }
        else if (matrix_sizes_[0] >= 3 && matrix_sizes_[1] >= 3
            && matrix_sizes_[2] >= 3)
        {
            return new EGFRDWorld(edge_lengths, matrix_sizes_);
        }
        else
        {
            return new EGFRDWorld(edge_lengths);
        }
    }

    virtual EGFRDWorld* create_world(const boost::shared_ptr<Model>& m) const
    {
        return extras::generate_world_from_model(*this, m);
    }

    virtual EGFRDSimulator* create_simulator(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<world_type>& world) const
    {
        if (user_max_shell_size_ > 0)
        {
            return new EGFRDSimulator(
                world, model, num_retries_, bd_dt_factor_, user_max_shell_size_);
        }
        else if (bd_dt_factor_ > 0)
        {
            return new EGFRDSimulator(
                world, model, num_retries_, bd_dt_factor_);
        }
        else if (num_retries_ >= 0)
        {
            return new EGFRDSimulator(world, model, num_retries_);
        }
        else
        {
            return new EGFRDSimulator(world, model);
        }
    }

    virtual EGFRDSimulator* create_simulator(
        const boost::shared_ptr<world_type>& world) const
    {
        if (boost::shared_ptr<Model> bound_model = world->lock_model())
        {
            return create_simulator(bound_model, world);
        }
        else
        {
            throw std::invalid_argument("A world must be bound to a model.");
        }
    }

protected:

    matrix_sizes_type matrix_sizes_;
    boost::shared_ptr<RandomNumberGenerator> rng_;

    Integer num_retries_;
    Integer bd_dt_factor_;
    Real user_max_shell_size_;
};

class BDFactory
    : public SimulatorFactory<EGFRDWorld, BDSimulator>
{
public:

    typedef SimulatorFactory<EGFRDWorld, BDSimulator> base_type;

protected:

    typedef EGFRDWorld::matrix_sizes_type matrix_sizes_type;

public:

    BDFactory(
        Integer dissociation_retry_moves=-1, Real bd_dt_factor=-1)
        : base_type(), matrix_sizes_(0, 0, 0), rng_(),
        num_retries_(dissociation_retry_moves),
        bd_dt_factor_(bd_dt_factor)
    {
        ; // do nothing
    }

    BDFactory(const matrix_sizes_type& matrix_sizes,
        Integer dissociation_retry_moves=-1, Real bd_dt_factor=-1)
        : base_type(), matrix_sizes_(matrix_sizes), rng_(),
        num_retries_(dissociation_retry_moves),
        bd_dt_factor_(bd_dt_factor)
    {
        ; // do nothing
    }

    BDFactory(const matrix_sizes_type& matrix_sizes,
        const boost::shared_ptr<RandomNumberGenerator>& rng,
        Integer dissociation_retry_moves=-1, Real bd_dt_factor=-1)
        : base_type(), matrix_sizes_(matrix_sizes), rng_(rng),
        num_retries_(dissociation_retry_moves),
        bd_dt_factor_(bd_dt_factor)
    {
        ; // do nothing
    }

    virtual ~BDFactory()
    {
        ; // do nothing
    }

    virtual EGFRDWorld* create_world(const std::string filename) const
    {
        return new EGFRDWorld(filename);
    }

    virtual EGFRDWorld* create_world(
        const Real3& edge_lengths = Real3(1, 1, 1)) const
    {
        if (rng_)
        {
            return new EGFRDWorld(edge_lengths, matrix_sizes_, rng_);
        }
        else if (matrix_sizes_[0] >= 3 && matrix_sizes_[1] >= 3
            && matrix_sizes_[2] >= 3)
        {
            return new EGFRDWorld(edge_lengths, matrix_sizes_);
        }
        else
        {
            return new EGFRDWorld(edge_lengths);
        }
    }

    virtual EGFRDWorld* create_world(const boost::shared_ptr<Model>& m) const
    {
        return extras::generate_world_from_model(*this, m);
    }

    virtual BDSimulator* create_simulator(
        const boost::shared_ptr<Model>& model,
        const boost::shared_ptr<world_type>& world) const
    {
        if (bd_dt_factor_ > 0)
        {
            return new BDSimulator(
                world, model, num_retries_, bd_dt_factor_);
        }
        else if (num_retries_ >= 0)
        {
            return new BDSimulator(world, model, num_retries_);
        }
        else
        {
            return new BDSimulator(world, model);
        }
    }

    virtual BDSimulator* create_simulator(
        const boost::shared_ptr<world_type>& world) const
    {
        if (boost::shared_ptr<Model> bound_model = world->lock_model())
        {
            return create_simulator(bound_model, world);
        }
        else
        {
            throw std::invalid_argument("A world must be bound to a model.");
        }
    }

protected:

    matrix_sizes_type matrix_sizes_;
    boost::shared_ptr<RandomNumberGenerator> rng_;

    Integer num_retries_;
    Integer bd_dt_factor_;
};

} // egfrd

} // ecell4

#endif /* __ECELL4_EGFRD_EGFRD_HPP */
