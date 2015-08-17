#ifndef __ECELL4_MESO_MESOSCOPIC_SIMULATOR_HPP
#define __ECELL4_MESO_MESOSCOPIC_SIMULATOR_HPP

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/SimulatorBase.hpp>
#include <ecell4/core/EventScheduler.hpp>

#include "MesoscopicWorld.hpp"

namespace ecell4
{

namespace meso
{

class MesoscopicSimulator
    : public SimulatorBase<Model, MesoscopicWorld>
{
public:

    typedef SimulatorBase<Model, MesoscopicWorld> base_type;
    typedef SubvolumeSpace::coordinate_type coordinate_type;

protected:

    class ReactionRuleProxyBase
    {
    public:

        ReactionRuleProxyBase()
            : sim_(), rr_()
        {
            ;
        }

        ReactionRuleProxyBase(MesoscopicSimulator* sim, const ReactionRule& rr)
            : sim_(sim), rr_(rr)
        {
            ;
        }

        virtual ~ReactionRuleProxyBase()
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
        virtual void inc(
            const Species& sp, const coordinate_type& c, const Integer val = +1) = 0;
        virtual const Real propensity(const coordinate_type& c) const = 0;

        inline void dec(const Species& sp, const coordinate_type& c)
        {
            inc(sp, c, -1);
        }

        virtual std::pair<ReactionRule, coordinate_type> draw(const coordinate_type& c) = 0;

    protected:

        inline const boost::shared_ptr<RandomNumberGenerator>& rng() const
        {
            return sim_->world()->rng();
        }

        inline const MesoscopicWorld& world() const
        {
            return (*sim_->world());
        }

    protected:

        MesoscopicSimulator* sim_;
        ReactionRule rr_;
    };

    class ReactionRuleProxy
        : public ReactionRuleProxyBase
    {
    public:

        typedef ReactionRuleProxyBase base_type;

        ReactionRuleProxy()
            : base_type()
        {
            ;
        }

        ReactionRuleProxy(MesoscopicSimulator* sim, const ReactionRule& rr)
            : base_type(sim, rr)
        {
            ;
        }

        virtual ~ReactionRuleProxy()
        {
            ;
        }

        std::pair<ReactionRule, coordinate_type> draw(const coordinate_type& c)
        {
            const std::pair<ReactionRule::reactant_container_type, Integer>
                retval(__draw(c));
            if (retval.second == 0)
            {
                return std::make_pair(ReactionRule(), c);
            }

            const std::vector<ReactionRule> reactions(generate(retval.first));

            assert(retval.second > 0);
            assert(retval.second >= reactions.size());

            if (reactions.size() == 0)
            {
                return std::make_pair(ReactionRule(), c);
            }
            else if (retval.second == 1)
            {
                // assert(possibles.size() == 1);
                return std::make_pair(reactions[0], c);
            }
            else
            {
                const std::vector<ReactionRule>::size_type
                    rnd2(static_cast<std::vector<ReactionRule>::size_type>(
                        rng()->uniform_int(0, retval.second - 1)));
                if (rnd2 >= reactions.size())
                {
                    return std::make_pair(ReactionRule(), c);
                }
                return std::make_pair(reactions[rnd2], c);
            }
        }

    protected:

        virtual std::pair<ReactionRule::reactant_container_type, Integer>
            __draw(const coordinate_type& c) = 0;
    };

    class DiffusionProxy
        : public ReactionRuleProxyBase
    {
    public:

        typedef ReactionRuleProxyBase base_type;

        DiffusionProxy()
            : base_type(), num_tot_()
        {
            ;
        }

        DiffusionProxy(MesoscopicSimulator* sim, const Species& sp)
            : base_type(sim, create_unimolecular_reaction_rule(sp, sp, 0.0)),
            num_tot_(sim->world()->num_subvolumes())
        {
            ;
        }

        virtual ~DiffusionProxy()
        {
            ;
        }

        std::pair<ReactionRule, coordinate_type> draw(const coordinate_type& c)
        {
            const Real3 lengths(sim_->world()->subvolume_edge_lengths());
            const Real px(1.0 / (lengths[0] * lengths[0])),
                py(1.0 / (lengths[1] * lengths[1])),
                pz(1.0 / (lengths[2] * lengths[2]));

            const Real rnd1(sim_->world()->rng()->uniform(0.0, px + py + pz));

            coordinate_type tgt;
            if (rnd1 < px * 0.5)
            {
                tgt = sim_->world()->get_neighbor(c, 0);
            }
            else if (rnd1 < px)
            {
                tgt = sim_->world()->get_neighbor(c, 1);
            }
            else if (rnd1 < px + py * 0.5)
            {
                tgt = sim_->world()->get_neighbor(c, 2);
            }
            else if (rnd1 < px + py)
            {
                tgt = sim_->world()->get_neighbor(c, 3);
            }
            else if (rnd1 < px + py + pz * 0.5)
            {
                tgt = sim_->world()->get_neighbor(c, 4);
            }
            else
            {
                tgt = sim_->world()->get_neighbor(c, 5);
            }
            return std::make_pair(rr_, tgt);
        }

        void initialize()
        {
            const Species sp(rr_.reactants()[0]);
            const Real D(sim_->world()->get_molecule_info(sp).D);
            const Real3 lengths(sim_->world()->subvolume_edge_lengths());
            const Real px(1.0 / (lengths[0] * lengths[0])),
                py(1.0 / (lengths[1] * lengths[1])),
                pz(1.0 / (lengths[2] * lengths[2]));
            k_ = 2 * D * (px + py + pz);

            std::fill(num_tot_.begin(), num_tot_.end(), 0);
            const std::vector<Species>& species(sim_->world()->species());
            for (std::vector<Species>::const_iterator i(species.begin());
                i != species.end(); ++i)
            {
                const Integer coef(get_coef(sp, *i));
                if (coef > 0)
                {
                    for (Integer j(0); j < sim_->world()->num_subvolumes(); ++j)
                    {
                        num_tot_[j] += coef * sim_->world()->num_molecules_exact(*i, j);
                    }
                }
            }
        }

        const Real propensity(const coordinate_type& c) const
        {
            return k_ * num_tot_[c];
        }

        void inc(const Species& sp, const coordinate_type& c, const Integer val = +1)
        {
            const Species& pttrn(rr_.reactants()[0]);
            const Integer coef(get_coef(pttrn, sp));
            if (coef > 0)
            {
                num_tot_[c] += coef * val;
            }
        }

    protected:

        std::vector<Real> num_tot_;
        Real k_;
    };

    class ZerothOrderReactionRuleProxy
        : public ReactionRuleProxy
    {
    public:

        typedef ReactionRuleProxy base_type;

        ZerothOrderReactionRuleProxy()
            : base_type()
        {
            ;
        }

        ZerothOrderReactionRuleProxy(MesoscopicSimulator* sim, const ReactionRule& rr)
            : base_type(sim, rr)
        {
            ;
        }

        void inc(const Species& sp, const coordinate_type& c, const Integer val = +1)
        {
            ; // do nothing
        }

        void initialize()
        {
            ; // do nothing
        }

        std::pair<ReactionRule::reactant_container_type, Integer> __draw(const coordinate_type& c)
        {
            return std::make_pair(ReactionRule::reactant_container_type(), 1);
        }

        const Real propensity(const coordinate_type& c) const
        {
            return rr_.k();
        }
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
            return (num_tot1_[c] * num_tot2_[c] - num_tot12_[c]) * rr_.k() / world().subvolume();
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
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());
            const ReactionRule::product_container_type& products(rr_.products());

            if (!(reactants.size() == 1 && products.size() == 1
                && reactants[0] == products[0] && coord_ == tgt_))
            {
                for (ReactionRule::product_container_type::const_iterator
                    it(rr_.products().begin());
                    it != rr_.products().end(); ++it)
                {
                    const Species& sp(*it);
                    if (!sim_->world()->on_structure(sp, tgt_))
                    {
                        // sim_->set_last_reaction(rr_);
                        update();
                        return; // XXX: do finalize
                    }
                }

                for (ReactionRule::reactant_container_type::const_iterator
                    it(rr_.reactants().begin());
                    it != rr_.reactants().end(); ++it)
                {
                    sim_->decrement_molecules(*it, coord_);
                }

                for (ReactionRule::product_container_type::const_iterator
                    it(rr_.products().begin());
                    it != rr_.products().end(); ++it)
                {
                    sim_->increment_molecules(*it, tgt_);
                }
            }

            if (coord_ != tgt_)
            {
                sim_->interrupt(tgt_);
                sim_->reset_last_reactions();
            }
            else
            {
                sim_->set_last_reaction(rr_);
            }

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
            dt_ = 0.0;
            while (true)
            {
                const std::pair<Real, std::pair<ReactionRule, coordinate_type> >
                    retval(sim_->draw_next_reaction(coord_));
                dt_ += retval.first;
                if (retval.first == inf || retval.second.first.k() > 0.0
                    || retval.second.second != coord_)
                {
                    rr_ = retval.second.first;
                    tgt_ = retval.second.second;
                    break;
                }
            }

            time_ += dt_;
        }

    protected:

        MesoscopicSimulator* sim_;
        coordinate_type coord_, tgt_;
        Real dt_;
        ReactionRule rr_;
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
    Real dt(void) const;
    Real next_time(void) const;

    void step(void);
    bool step(const Real & upto);

    // Optional members
    std::vector<ReactionRule> last_reactions() const;
    void set_last_reaction(const ReactionRule& rr);

    void reset_last_reactions()
    {
        last_reactions_.clear();
    }

    /**
     * recalculate reaction propensities and draw the next time.
     */
    void initialize();

    inline boost::shared_ptr<RandomNumberGenerator> rng()
    {
        return (*world_).rng();
    }

    void interrupt(const coordinate_type& coord)
    {
        interrupted_ = coord;
    }

protected:

    void interrupt_all(const Real& t);
    std::pair<Real, std::pair<ReactionRule, coordinate_type> >
        draw_next_reaction(const coordinate_type& c);
    std::pair<Real, std::pair<ReactionRule, coordinate_type> >
        draw_next_diffusion(const coordinate_type& c);
    void increment_molecules(const Species& sp, const coordinate_type& c);
    void decrement_molecules(const Species& sp, const coordinate_type& c);

protected:

    std::vector<ReactionRule> last_reactions_;

    boost::ptr_vector<ReactionRuleProxyBase> proxies_;
    EventScheduler scheduler_;
    std::vector<EventScheduler::identifier_type> event_ids_;
    coordinate_type interrupted_;
};

} // meso

} // ecell4

#endif /* __ECELL4_MESO_MESOSCOPIC_SIMULATOR_HPP */
