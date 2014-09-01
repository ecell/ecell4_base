#ifndef __ECELL4_MESO_MESOSCOPIC_SIMULATOR_HPP
#define __ECELL4_MESO_MESOSCOPIC_SIMULATOR_HPP

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/Simulator.hpp>
#include <ecell4/core/EventScheduler.hpp>

#include "MesoscopicWorld.hpp"

namespace ecell4
{

namespace meso
{

class MesoscopicSimulator
    : public Simulator<Model, MesoscopicWorld>
{
public:

    typedef Simulator<Model, MesoscopicWorld> base_type;
    typedef SubvolumeSpace::coordinate_type coordinate_type;

protected:

    class ReactionRuleProxy
    {
    public:

        ReactionRuleProxy()
            : sim_(), rr_()
        {
            ;
        }

        ReactionRuleProxy(MesoscopicSimulator* sim, const ReactionRule& rr)
            : sim_(sim), rr_(rr)
        {
            ;
        }

        virtual ~ReactionRuleProxy()
        {
            ;
        }

        inline const Integer get_coef(const Species& pttrn, const Species& sp) const
        {
            return sim_->model()->apply(pttrn, sp);
        }

        inline const std::vector<ReactionRule> generate(
            const ReactionRule::reactant_container_type& reactants) const
        {
            return sim_->model()->apply(rr_, reactants);
        }

        virtual void initialize() = 0;
        virtual void inc(const Species& sp, const coordinate_type& c, const Integer val = +1) = 0;
        virtual const Real propensity(const coordinate_type& c) const = 0;

        inline void dec(const Species& sp, const coordinate_type& c)
        {
            inc(sp, c, -1);
        }

        ReactionRule draw(const coordinate_type& c)
        {
            const std::pair<ReactionRule::reactant_container_type, Integer>
                retval(__draw(c));
            if (retval.second == 0)
            {
                return ReactionRule();
            }

            const std::vector<ReactionRule> reactions(generate(retval.first));

            assert(retval.second > 0);
            assert(retval.second >= reactions.size());

            if (reactions.size() == 0)
            {
                return ReactionRule();
            }
            else if (retval.second == 1)
            {
                // assert(possibles.size() == 1);
                return reactions[0];
            }
            else
            {
                const Integer rnd2(rng()->uniform_int(0, retval.second - 1));
                if (rnd2 >= reactions.size())
                {
                    return ReactionRule();
                }
                return reactions[rnd2];
            }
        }

    protected:

        inline const boost::shared_ptr<RandomNumberGenerator>& rng() const
        {
            return sim_->world()->rng();
        }

        inline const MesoscopicWorld& world() const
        {
            return (*sim_->world());
        }

        virtual std::pair<ReactionRule::reactant_container_type, Integer>
            __draw(const coordinate_type& c) = 0;

    protected:

        MesoscopicSimulator* sim_;
        ReactionRule rr_;
    };

    class FirstOrderReactionRuleProxy
        : public ReactionRuleProxy
    {
    public:

        typedef ReactionRuleProxy base_type;

        FirstOrderReactionRuleProxy()
            : base_type(), num_tot1_()
        {
            ;
        }

        FirstOrderReactionRuleProxy(MesoscopicSimulator* sim, const ReactionRule& rr)
            : base_type(sim, rr), num_tot1_(sim->world()->num_subvolumes())
        {
            ;
        }

        void inc(const Species& sp, const coordinate_type& c, const Integer val = +1)
        {
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());
            const Integer coef(get_coef(reactants[0], sp));
            if (coef > 0)
            {
                num_tot1_[c] += coef * val;
            }
        }

        void initialize()
        {
            const std::vector<Species>& species(world().list_species());
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());

            std::fill(num_tot1_.begin(), num_tot1_.end(), 0);
            for (std::vector<Species>::const_iterator i(species.begin());
                i != species.end(); ++i)
            {
                const Integer coef(get_coef(reactants[0], *i));
                if (coef > 0)
                {
                    for (Integer j(0); j < world().num_subvolumes(); ++j)
                    {
                        num_tot1_[j] += coef * world().num_molecules_exact(*i, j);
                    }
                }
            }
        }

        std::pair<ReactionRule::reactant_container_type, Integer> __draw(const coordinate_type& c)
        {
            const std::vector<Species>& species(world().list_species());
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());

            const Real rnd1(rng()->uniform(0.0, num_tot1_[c]));

            Integer num_tot(0);
            for (std::vector<Species>::const_iterator i(species.begin());
                i != species.end(); ++i)
            {
                const Integer coef(get_coef(reactants[0], *i));
                if (coef > 0)
                {
                    num_tot += coef * world().num_molecules_exact(*i, c);
                    if (num_tot >= rnd1)
                    {
                        return std::make_pair(
                            ReactionRule::reactant_container_type(1, *i), coef);
                    }
                }
            }

            return std::make_pair(ReactionRule::reactant_container_type(), 0);
        }

        const Real propensity(const coordinate_type& c) const
        {
            return num_tot1_[c] * rr_.k();
        }

    protected:

        std::vector<Integer> num_tot1_;
    };

    class SecondOrderReactionRuleProxy:
        public ReactionRuleProxy
    {
    public:

        typedef ReactionRuleProxy base_type;

        SecondOrderReactionRuleProxy()
            : base_type(), num_tot1_(0), num_tot2_(0), num_tot12_(0)
        {
            ;
        }

        SecondOrderReactionRuleProxy(MesoscopicSimulator* sim, const ReactionRule& rr)
            : base_type(sim, rr), num_tot1_(sim->world()->num_subvolumes()), num_tot2_(sim->world()->num_subvolumes()), num_tot12_(sim->world()->num_subvolumes())
        {
            ;
        }

        void inc(const Species& sp, const coordinate_type& c, const Integer val = +1)
        {
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());
            const Integer coef1(get_coef(reactants[0], sp));
            const Integer coef2(get_coef(reactants[1], sp));
            if (coef1 > 0 || coef2 > 0)
            {
                const Integer tmp(coef1 * val);
                num_tot1_[c] += tmp;
                num_tot2_[c] += coef2 * val;
                num_tot12_[c] += coef2 * tmp;
            }
        }

        void initialize()
        {
            const std::vector<Species>& species(world().list_species());
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());

            std::fill(num_tot1_.begin(), num_tot1_.end(), 0);
            std::fill(num_tot2_.begin(), num_tot2_.end(), 0);
            std::fill(num_tot12_.begin(), num_tot12_.end(), 0);
            for (std::vector<Species>::const_iterator i(species.begin());
                i != species.end(); ++i)
            {
                const Integer coef1(get_coef(reactants[0], *i));
                const Integer coef2(get_coef(reactants[1], *i));
                if (coef1 > 0 || coef2 > 0)
                {
                    for (Integer j(0); j < world().num_subvolumes(); ++j)
                    {
                        const Integer num(world().num_molecules_exact(*i, j));
                        const Integer tmp(coef1 * num);
                        num_tot1_[j] += tmp;
                        num_tot2_[j] += coef2 * num;
                        num_tot12_[j] += coef2 * tmp;
                    }
                }
            }
        }

        std::pair<ReactionRule::reactant_container_type, Integer> __draw(const coordinate_type& c)
        {
            const std::vector<Species>& species(world().list_species());
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());

            const Real rnd1(rng()->uniform(0.0, num_tot1_[c]));

            Integer num_tot(0), coef1(0);
            std::vector<Species>::const_iterator itr1(species.begin());
            for (; itr1 != species.end(); ++itr1)
            {
                const Integer coef(get_coef(reactants[0], *itr1));
                if (coef > 0)
                {
                    num_tot += coef * world().num_molecules_exact(*itr1, c);
                    if (num_tot >= rnd1)
                    {
                        coef1 = coef;
                        break;
                    }
                }
            }

            const Real rnd2(
                rng()->uniform(0.0, num_tot2_[c] - get_coef(reactants[0], *itr1)));

            num_tot = 0;
            for (std::vector<Species>::const_iterator i(species.begin());
                i != species.end(); ++i)
            {
                const Integer coef(get_coef(reactants[1], *i));
                if (coef > 0)
                {
                    const Integer num(world().num_molecules_exact(*i, c));
                    num_tot += coef * (i == itr1 ? num - 1 : num);
                    if (num_tot >= rnd2)
                    {
                        ReactionRule::reactant_container_type exact_reactants(2);
                        exact_reactants[0] = *itr1;
                        exact_reactants[1] = *i;
                        return std::make_pair(exact_reactants, coef1 * coef);
                    }
                }
            }

            return std::make_pair(ReactionRule::reactant_container_type(), 0);
        }

        const Real propensity(const coordinate_type& c) const
        {
            return (num_tot1_[c] * num_tot2_[c] - num_tot12_[c]) * rr_.k() / world().volume();
        }

    protected:

        std::vector<Integer> num_tot1_, num_tot2_, num_tot12_;
    };

    struct SubvolumeEvent
        : public EventScheduler::Event
    {
    public:

        SubvolumeEvent(MesoscopicSimulator* sim, const coordinate_type& c, const Real& t)
            : EventScheduler::Event(t), sim_(sim), coord_(c)
        {
            update();
        }

        virtual ~SubvolumeEvent()
        {
            ;
        }

        virtual void fire()
        {
            for (ReactionRule::reactant_container_type::const_iterator
                it(next_.reactants().begin());
                it != next_.reactants().end(); ++it)
            {
                sim_->decrement_molecules(*it, coord_);
            }

            for (ReactionRule::product_container_type::const_iterator
                it(next_.products().begin());
                it != next_.products().end(); ++it)
            {
                sim_->increment_molecules(*it, coord_);
            }

            sim_->set_last_reaction(next_);

            update();
        }

        virtual void interrupt(Real const& t)
        {
            time_ = t;
            update();
        }

        coordinate_type coordinate() const
        {
            return coord_;
        }

        void update()
        {
            while (true)
            {
                const std::pair<Real, ReactionRule>
                    retval(sim_->draw_next_reaction(coord_));
                time_ += retval.first;
                if (retval.first == inf || retval.second.k() > 0.0)
                {
                    next_ = retval.second;
                    break;
                }
            }
        }

    protected:

        MesoscopicSimulator* sim_;
        coordinate_type coord_;
        ReactionRule next_;
    };

public:

    MesoscopicSimulator(
        boost::shared_ptr<Model> model,
        boost::shared_ptr<MesoscopicWorld> world)
        : base_type(model, world)
    {
        initialize();
    }

    MesoscopicSimulator(boost::shared_ptr<MesoscopicWorld> world)
        : base_type(world)
    {
        initialize();
    }

    // SimulatorTraits

    Real t(void) const;
    Real dt(void) const;
    Real next_time(void) const;

    void step(void) ;
    bool step(const Real & upto);

    // Optional members

    void set_t(const Real &t);
    std::vector<ReactionRule> last_reactions() const;
    void set_last_reaction(const ReactionRule& rr);

    /**
     * recalculate reaction propensities and draw the next time.
     */
    void initialize();

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return (*world_).rng();
    }

protected:

    void interrupt(const Real& t);
    std::pair<Real, ReactionRule> draw_next_reaction(const coordinate_type& c);
    void increment_molecules(const Species& sp, const coordinate_type& c);
    void decrement_molecules(const Species& sp, const coordinate_type& c);

protected:

    std::vector<ReactionRule> last_reactions_;

    boost::ptr_vector<ReactionRuleProxy> proxies_;
    EventScheduler scheduler_;
};

} // meso

} // ecell4

#endif /* __ECELL4_MESO_MESOSCOPIC_SIMULATOR_HPP */
