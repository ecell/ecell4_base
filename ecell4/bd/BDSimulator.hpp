#ifndef ECELL4_BD_BD_SIMULATOR_HPP
#define ECELL4_BD_BD_SIMULATOR_HPP

#include <stdexcept>
#include <boost/shared_ptr.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/SimulatorBase.hpp>

#include "BDWorld.hpp"
#include "BDPropagator.hpp"
#include "BDPropagator2D.hpp"


namespace ecell4
{

namespace bd
{

class BDSimulator
    : public SimulatorBase<Model, BDWorld>
{
public:

    typedef SimulatorBase<Model, BDWorld> base_type;
    typedef BDPropagator::reaction_info_type reaction_info_type;

public:

    BDSimulator(boost::shared_ptr<Model> model,
        boost::shared_ptr<BDWorld> world, Real bd_dt_factor = 1e-5)
        : base_type(model, world), dt_(0), bd_dt_factor_(bd_dt_factor)
    {
        initialize();
    }

    BDSimulator(boost::shared_ptr<BDWorld> world, Real bd_dt_factor = 1e-5)
        : base_type(world), dt_(0), bd_dt_factor_(bd_dt_factor)
    {
        initialize();
    }

    // SimulatorTraits

    void initialize()
    {
        last_reactions_.clear();
        dt_ = determine_dt();

        // XXX: refine determine_reaction_length!
        reaction_length_ = this->min_sigma() * 1e-2;
    }

    Real determine_dt() const
    {
        const std::vector<Species> splist(world_->list_species());

        Real rmin(inf), Dmax(0.0);
        for (std::vector<Species>::const_iterator i(splist.begin());
            i != splist.end(); ++i)
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
    }

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return (*world_).rng();
    }
private:

    Real min_sigma()
    {
        std::vector<Species> sps(world_->list_species());
        Real minimum(std::numeric_limits<Real>::max());
        for(std::vector<Species>::const_iterator
                iter = sps.begin(); iter != sps.end(); ++iter)
        {
            const BDWorld::molecule_info_type mol(
                    world_->get_molecule_info(*iter));
            if(minimum > mol.radius) minimum = mol.radius;
        }
        return minimum;
    }

protected:

    /**
     * the protected internal state of BDSimulator.
     * they are needed to be saved/loaded with Visitor pattern.
     */
    Real dt_;
    Real reaction_length_;
    const Real bd_dt_factor_;
    std::vector<std::pair<ReactionRule, reaction_info_type> > last_reactions_;
};

} // bd

} // ecell4

#endif /* ECELL4_BD_BD_SIMULATOR_HPP */
