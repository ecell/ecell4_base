#ifndef ECELL4_EGFRD_EGFRD_HPP
#define ECELL4_EGFRD_EGFRD_HPP

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

typedef EGFRDSimulator::reaction_info_type ReactionInfo;
// typedef BDSimulator::reaction_info_type ReactionInfo;

class EGFRDFactory
    : public SimulatorFactory<EGFRDWorld, EGFRDSimulator>
{
public:

    typedef SimulatorFactory<EGFRDWorld, EGFRDSimulator> base_type;

protected:

    typedef EGFRDWorld::matrix_sizes_type matrix_sizes_type;

public:

    EGFRDFactory(
        const matrix_sizes_type& matrix_sizes = default_matrix_sizes(),
        Real bd_dt_factor = default_bd_dt_factor(),
        Integer dissociation_retry_moves = default_dissociation_retry_moves(),
        Real user_max_shell_size = default_user_max_shell_size())
        : base_type(), rng_(),
          matrix_sizes_(matrix_sizes), bd_dt_factor_(bd_dt_factor),
          dissociation_retry_moves_(dissociation_retry_moves),
          user_max_shell_size_(user_max_shell_size)
    {
        ; // do nothing
    }

    static inline const matrix_sizes_type default_matrix_sizes()
    {
        return Integer3(0, 0, 0);
    }

    static inline const Real default_bd_dt_factor()
    {
        return 0.0;
    }

    static inline const Integer default_dissociation_retry_moves()
    {
        return -1;
    }

    static inline const Real default_user_max_shell_size()
    {
        return 0.0;
    }

    virtual ~EGFRDFactory()
    {
        ; // do nothing
    }

    EGFRDFactory& rng(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        rng_ = rng;
        return (*this);
    }

    inline EGFRDFactory* rng_ptr(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        return &(this->rng(rng));  //XXX: == this
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
        // else if (matrix_sizes_[0] >= 3 && matrix_sizes_[1] >= 3
        //     && matrix_sizes_[2] >= 3)
        else if (matrix_sizes_ != default_matrix_sizes())
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
        if (user_max_shell_size_ != default_user_max_shell_size())
        {
            return new EGFRDSimulator(
                world, model, bd_dt_factor_, dissociation_retry_moves_, user_max_shell_size_);
        }
        else if (dissociation_retry_moves_ != default_dissociation_retry_moves())
        {
            return new EGFRDSimulator(
                world, model, bd_dt_factor_, dissociation_retry_moves_);
        }
        else if (bd_dt_factor_ != default_bd_dt_factor())
        {
            return new EGFRDSimulator(world, model, bd_dt_factor_);
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

    boost::shared_ptr<RandomNumberGenerator> rng_;
    matrix_sizes_type matrix_sizes_;
    Real bd_dt_factor_;
    Integer dissociation_retry_moves_;
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
        const matrix_sizes_type& matrix_sizes = default_matrix_sizes(),
        Real bd_dt_factor = default_bd_dt_factor(),
        Integer dissociation_retry_moves = default_dissociation_retry_moves())
        : base_type(), rng_(),
          matrix_sizes_(matrix_sizes), bd_dt_factor_(bd_dt_factor),
          dissociation_retry_moves_(dissociation_retry_moves)
    {
        ; // do nothing
    }

    static inline const matrix_sizes_type default_matrix_sizes()
    {
        return Integer3(0, 0, 0);
    }

    static inline const Real default_bd_dt_factor()
    {
        return 0.0;
    }

    static inline const Integer default_dissociation_retry_moves()
    {
        return -1;
    }

    virtual ~BDFactory()
    {
        ; // do nothing
    }

    BDFactory& rng(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        rng_ = rng;
        return (*this);
    }

    inline BDFactory* rng_ptr(const boost::shared_ptr<RandomNumberGenerator>& rng)
    {
        return &(this->rng(rng));  //XXX: == this
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
        // else if (matrix_sizes_[0] >= 3 && matrix_sizes_[1] >= 3
        //     && matrix_sizes_[2] >= 3)
        else if (matrix_sizes_ != default_matrix_sizes())
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
        if (dissociation_retry_moves_ != default_dissociation_retry_moves())
        {
            return new BDSimulator(
                world, model, bd_dt_factor_, dissociation_retry_moves_);
        }
        else if (bd_dt_factor_ != default_bd_dt_factor())
        {
            return new BDSimulator(world, model, bd_dt_factor_);
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

    boost::shared_ptr<RandomNumberGenerator> rng_;
    matrix_sizes_type matrix_sizes_;
    Real bd_dt_factor_;
    Integer dissociation_retry_moves_;
};

} // egfrd

} // ecell4

#endif /* ECELL4_EGFRD_EGFRD_HPP */
