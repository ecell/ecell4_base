#ifndef ECELL4_BD_BD_SIMULATOR_HPP
#define ECELL4_BD_BD_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/SimulatorBase.hpp>

#include "BDWorld.hpp"
#include "BDPropagator.hpp"


namespace ecell4
{

namespace bd
{

class BDSimulator
    : public SimulatorBase<BDWorld>
{
public:

    typedef SimulatorBase<BDWorld> base_type;
    typedef BDPropagator::reaction_info_type reaction_info_type;

public:

    BDSimulator(
        boost::shared_ptr<BDWorld> world, boost::shared_ptr<Model> model,
        Real bd_dt_factor = 1e-5)
        : base_type(world, model), dt_(0), bd_dt_factor_(bd_dt_factor), dt_set_by_user_(false)
    {
        initialize();
    }

    BDSimulator(boost::shared_ptr<BDWorld> world, Real bd_dt_factor = 1e-5)
        : base_type(world), dt_(0), bd_dt_factor_(bd_dt_factor), dt_set_by_user_(false)
    {
        initialize();
    }

    // SimulatorTraits

    void initialize()
    {
        last_reactions_.clear();
        if (!dt_set_by_user_)
        {
            dt_ = determine_dt();
        }
    }

    Real determine_dt() const
    {
        constexpr Real inf = std::numeric_limits<Real>::infinity();
        Real rmin(inf), Dmax(0.0);

        for (std::vector<Species>::const_iterator i(model_->species_attributes().begin());
            i != model_->species_attributes().end(); ++i)
        {
            const BDWorld::molecule_info_type
                info(world_->get_molecule_info(*i));

            if (rmin > info.radius)
            {
                rmin = info.radius;
            }
            if (Dmax < info.D)
            {
                Dmax = info.D;
            }
        }

        // const std::vector<Species> splist(world_->list_species());

        // for (std::vector<Species>::const_iterator i(splist.begin());
        //     i != splist.end(); ++i)
        // {
        //     const BDWorld::molecule_info_type
        //         info(world_->get_molecule_info(*i));
        //     if (rmin > info.radius)
        //     {
        //         rmin = info.radius;
        //     }
        //     if (Dmax < info.D)
        //     {
        //         Dmax = info.D;
        //     }
        // }

        const Real dt(rmin < inf && Dmax > 0.0
            ? 4.0 * rmin * rmin / (2.0 * Dmax) * bd_dt_factor_
            // ? rmin * rmin / (6.0 * Dmax) * bd_dt_factor_
            : inf);
        return dt;
    }

    Real dt() const
    {
        return dt_;
    }

    void step();
    bool step(const Real& upto);

    // Optional members

    virtual bool check_reaction() const
    {
        return last_reactions_.size() > 0;
    }

    std::vector<std::pair<ReactionRule, reaction_info_type> >
        last_reactions() const
    {
        return last_reactions_;
    }

    void set_dt(const Real& dt)
    {
        if (dt <= 0)
        {
            throw std::invalid_argument("The step size must be positive.");
        }
        dt_ = dt;
        dt_set_by_user_ = true;
    }

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return (*world_).rng();
    }

protected:

    void attempt_synthetic_reaction(const ReactionRule& rr);

protected:

    /**
     * the protected internal state of BDSimulator.
     * they are needed to be saved/loaded with Visitor pattern.
     */
    Real dt_;
    const Real bd_dt_factor_;
    bool dt_set_by_user_;
    std::vector<std::pair<ReactionRule, reaction_info_type> > last_reactions_;
};

} // bd

} // ecell4

#endif /* ECELL4_BD_BD_SIMULATOR_HPP */
