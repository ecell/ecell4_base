#ifndef ECELL4_MESO_MESOSCOPIC_SIMULATOR_HPP
#define ECELL4_MESO_MESOSCOPIC_SIMULATOR_HPP

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

class ReactionInfo
{
public:

    typedef SubvolumeSpace::coordinate_type coordinate_type;
    typedef Species element_type;
    typedef std::vector<element_type> container_type;

public:

    ReactionInfo(
        const Real t,
        const container_type& reactants,
        const container_type& products,
        const coordinate_type coord)
        : t_(t), reactants_(reactants), products_(products), coord_(coord)
    {}

    ReactionInfo(
        const ReactionInfo& another)
        : t_(another.t()), reactants_(another.reactants()), products_(another.products()),
          coord_(another.coordinate())
    {}

    Real t() const
    {
        return t_;
    }

    const container_type& reactants() const
    {
        return reactants_;
    }

    void add_reactant(const element_type& elem)
    {
        reactants_.push_back(elem);
    }

    const container_type& products() const
    {
        return products_;
    }

    void add_product(const element_type& elem)
    {
        products_.push_back(elem);
    }

    coordinate_type coordinate() const
    {
        return coord_;
    }

protected:

    Real t_;
    container_type reactants_, products_;
    coordinate_type coord_;
};

class MesoscopicSimulator
    : public SimulatorBase<Model, MesoscopicWorld>
{
public:

    typedef SimulatorBase<Model, MesoscopicWorld> base_type;
    typedef SubvolumeSpace::coordinate_type coordinate_type;
    typedef ReactionInfo reaction_info_type;

protected:

    class ReactionRuleProxyBase
    {
    public:

        ReactionRuleProxyBase()
            : sim_()
        {
            ;
        }

        ReactionRuleProxyBase(MesoscopicSimulator* sim)
            : sim_(sim)
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

        virtual void initialize() = 0;
        virtual const Real propensity(const coordinate_type& c) const = 0;
        virtual void fire(const Real t, const coordinate_type& src) = 0;

        virtual void inc(
            const Species& sp, const coordinate_type& c, const Integer val = +1) = 0;

        inline void dec(const Species& sp, const coordinate_type& c)
        {
            inc(sp, c, -1);
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

    protected:

        MesoscopicSimulator* sim_;
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
            : base_type(sim), rr_(rr)
        {
            ;
        }

        virtual ~ReactionRuleProxy()
        {
            ;
        }

        const ReactionRule& reaction_rule() const
        {
            return rr_;
        }

        std::vector<Integer> check_dependency(const Species& sp) const
        {
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());
            std::vector<Integer> coefs(reactants.size(), 0);
            for (std::vector<Integer>::size_type i = 0; i < coefs.size(); ++i)
            {
                coefs[i] = get_coef(reactants[i], sp);
            }
            return coefs;
        }

        virtual void inc_with_coefs(const std::vector<Integer>& coefs,
                         const coordinate_type& c, const Integer val = +1) = 0;

        inline const std::vector<ReactionRule> generate(
            const ReactionRule::reactant_container_type& reactants) const
        {
            return sim_->model()->apply(rr_, reactants);
        }

        std::pair<ReactionRule, coordinate_type> draw(const coordinate_type& c)
        {
            const std::pair<ReactionRule::reactant_container_type, Integer>
                retval(__draw(c));
            const std::vector<ReactionRule> reactions(generate(retval.first));

            assert(retval.second > 0);
            assert(retval.second >= reactions.size());

            if (reactions.size() == 0)
            {
                return std::make_pair(ReactionRule(), c);
            }
            else if (retval.second == 1)
            {
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

        virtual void fire(const Real t, const coordinate_type& src)
        {
            const std::pair<ReactionRule, coordinate_type>
                retval = this->draw(src);

            const ReactionRule& nextr = retval.first;
            // const coordinate_type& dst = retval.second;
            const ReactionRule::reactant_container_type& reactants(nextr.reactants());
            const ReactionRule::product_container_type& products(nextr.products());

            assert(retval.second == src);

            for (ReactionRule::product_container_type::const_iterator
                    it(products.begin()); it != products.end(); ++it)
            {
                const Species& sp(*it);

                if (!sim_->world()->on_structure(sp, src))
                {
                    ; // do nothing except for update()
                    return;
                }
            }

            for (ReactionRule::reactant_container_type::const_iterator
                    it(reactants.begin()); it != reactants.end(); ++it)
            {
                sim_->decrement_molecules(*it, src);
            }

            for (ReactionRule::product_container_type::const_iterator
                    it(products.begin()); it != products.end(); ++it)
            {
                sim_->increment_molecules(*it, src);
            }

            sim_->add_last_reaction(
                nextr, reaction_info_type(t, reactants, products, src));
        }

    protected:

        virtual std::pair<ReactionRule::reactant_container_type, Integer>
            __draw(const coordinate_type& c) = 0;

    protected:

        ReactionRule rr_;
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

        void inc_with_coefs(const std::vector<Integer>& coefs,
                 const coordinate_type& c, const Integer val = +1)
        {
            ; // do nothing
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
            return rr_.k() * world().subvolume();
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

        void inc_with_coefs(const std::vector<Integer>& coefs,
                 const coordinate_type& c, const Integer val = +1)
        {
            num_tot1_[c] += coefs[0] * val;
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

            throw IllegalState("FirstOrderReactionRuleProxy: Never reach here.");
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

        void inc_with_coefs(const std::vector<Integer>& coefs,
                 const coordinate_type& c, const Integer val = +1)
        {
            num_tot1_[c] += coefs[0] * val;
            num_tot2_[c] += coefs[1] * val;
            num_tot12_[c] += coefs[0] * coefs[1] * val;
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

            throw IllegalState("SecondOrderReactionRuleProxy: Never reach here.");
        }

        const Real propensity(const coordinate_type& c) const
        {
            return (num_tot1_[c] * num_tot2_[c] - num_tot12_[c]) * rr_.k() / world().subvolume();
        }

    protected:

        std::vector<Integer> num_tot1_, num_tot2_, num_tot12_;
    };

    class StructureSecondOrderReactionRuleProxy:
        public ReactionRuleProxy
    {
    public:

        typedef ReactionRuleProxy base_type;

        StructureSecondOrderReactionRuleProxy()
            : base_type(), num_tot_(0), stidx_(0), spidx_(0)
        {
            ;
        }

        StructureSecondOrderReactionRuleProxy(
            MesoscopicSimulator* sim, const ReactionRule& rr,
            const ReactionRule::reactant_container_type::size_type stidx)
            : base_type(sim, rr), num_tot_(sim->world()->num_subvolumes()),
              stidx_(stidx), spidx_(stidx == 0 ? 1 : 0)
        {
            ;
        }

        void inc_with_coefs(const std::vector<Integer>& coefs,
                 const coordinate_type& c, const Integer val = +1)
        {
            num_tot_[c] += coefs[spidx_] * val;
        }

        void inc(const Species& sp, const coordinate_type& c, const Integer val = +1)
        {
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());
            const Integer coef(get_coef(reactants[spidx_], sp));
            if (coef > 0)
            {
                num_tot_[c] += coef * val;
            }
        }

        void initialize()
        {
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());

            // if (world().get_dimension(reactants[stidx_]) != Shape::TWO)
            // {
            //     throw NotSupported(
            //         "A second order reaction is only acceptable"
            //         " with a structure with dimension two.");
            // } else
            if (world().has_structure(reactants[spidx_]))
            {
                throw NotSupported(
                    "A second order reaction between structures has no mean.");
            }

            const std::vector<Species>& species(world().list_species());
            std::fill(num_tot_.begin(), num_tot_.end(), 0);
            for (std::vector<Species>::const_iterator i(species.begin());
                i != species.end(); ++i)
            {
                const Integer coef(get_coef(reactants[spidx_], *i));
                if (coef > 0)
                {
                    for (Integer j(0); j < world().num_subvolumes(); ++j)
                    {
                        num_tot_[j] += coef * world().num_molecules_exact(*i, j);
                    }
                }
            }
        }

        std::pair<ReactionRule::reactant_container_type, Integer> __draw(const coordinate_type& c)
        {
            const std::vector<Species>& species(world().list_species());
            const ReactionRule::reactant_container_type& reactants(rr_.reactants());

            const Real rnd1(rng()->uniform(0.0, num_tot_[c]));

            Integer tot(0);
            for (std::vector<Species>::const_iterator i(species.begin());
                i != species.end(); ++i)
            {
                const Integer coef(get_coef(reactants[spidx_], *i));
                if (coef > 0)
                {
                    tot += coef * world().num_molecules_exact(*i, c);
                    if (tot >= rnd1)
                    {
                        ReactionRule::reactant_container_type retval(2);
                        retval[spidx_] = *i;
                        retval[stidx_] = reactants[stidx_];
                        return std::make_pair(retval, coef);
                    }
                }
            }
            throw IllegalState("StructureSecondOrderReactionRuleProxy: Never reach here.");
        }

        const Real propensity(const coordinate_type& c) const
        {
            return (num_tot_[c] * rr_.k()
                    * world().get_occupancy(rr_.reactants()[stidx_], c));
        }

    protected:

        std::vector<Integer> num_tot_;
        ReactionRule::reactant_container_type::size_type stidx_, spidx_;
    };

    class DiffusionProxy
        : public ReactionRuleProxyBase
    {
    public:

        typedef ReactionRuleProxyBase base_type;

    protected:

        typedef std::pair<ReactionRuleProxy*, std::vector<Integer> >
            dependency_type;
        typedef std::vector<dependency_type> dependency_container_type;

    public:

        DiffusionProxy()
            : base_type(), pool_(), dependencies_()
        {
            ;
        }

        DiffusionProxy(MesoscopicSimulator* sim, const Species& sp)
            : base_type(sim), pool_(sim->world()->get_pool(sp)), dependencies_()
        {
            ;
        }

        virtual ~DiffusionProxy()
        {
            ;
        }

        coordinate_type draw(const coordinate_type& c)
        {
            const Real3 lengths(sim_->world()->subvolume_edge_lengths());
            const Real px(1.0 / (lengths[0] * lengths[0])),
                py(1.0 / (lengths[1] * lengths[1])),
                pz(1.0 / (lengths[2] * lengths[2]));

            const Real rnd1(sim_->world()->rng()->uniform(0.0, px + py + pz));

            if (rnd1 < px * 0.5)
            {
                return sim_->world()->get_neighbor(c, 0);
            }
            else if (rnd1 < px)
            {
                return sim_->world()->get_neighbor(c, 1);
            }
            else if (rnd1 < px + py * 0.5)
            {
                return sim_->world()->get_neighbor(c, 2);
            }
            else if (rnd1 < px + py)
            {
                return sim_->world()->get_neighbor(c, 3);
            }
            else if (rnd1 < px + py + pz * 0.5)
            {
                return sim_->world()->get_neighbor(c, 4);
            }
            return sim_->world()->get_neighbor(c, 5);
        }

        void initialize()
        {
            const Real D = pool_->D();
            const Real3 lengths(sim_->world()->subvolume_edge_lengths());
            const Real px(1.0 / (lengths[0] * lengths[0])),
                py(1.0 / (lengths[1] * lengths[1])),
                pz(1.0 / (lengths[2] * lengths[2]));
            k_ = 2 * D * (px + py + pz);
        }

        const Real propensity(const coordinate_type& c) const
        {
            return k_ * pool_->num_molecules(c);
        }

        void inc(const Species& sp, const coordinate_type& c, const Integer val = +1)
        {
            ; // do nothing
        }

        virtual void fire(const Real t, const coordinate_type& src)
        {
            const coordinate_type dst = this->draw(src);

            if (dst == src)
            {
                ; // do nothing except for update()
                return;
            }

            if (!sim_->world()->check_structure(pool_->loc(), dst))
            {
                ; // do nothing except for update()
                return;
            }

            {
                // sim_->decrement(pool_, src);
                // sim_->increment(pool_, dst);

                pool_->remove_molecules(1, src);
                pool_->add_molecules(1, dst);

                for (dependency_container_type::const_iterator i(dependencies_.begin());
                     i != dependencies_.end(); ++i)
                {
                    (*i).first->inc_with_coefs((*i).second, src, -1);
                    (*i).first->inc_with_coefs((*i).second, dst, +1);
                }
            }

            sim_->interrupt(dst);
        }

        void set_dependency(ReactionRuleProxy* proxy)
        {
            const std::vector<Integer> coefs = proxy->check_dependency(pool_->species());
            if (std::count(coefs.begin(), coefs.end(), 0) < coefs.size())
            {
                dependencies_.push_back(std::make_pair(proxy, coefs));
            }
        }

    protected:

        const boost::shared_ptr<MesoscopicWorld::PoolBase> pool_;
        Real k_;

        dependency_container_type dependencies_;
    };

    struct SubvolumeEvent
        : public Event
    {
    public:

        SubvolumeEvent(MesoscopicSimulator* sim, const coordinate_type& c, const Real& t)
            : Event(t), sim_(sim), coord_(c)
        {
            update();
        }

        virtual ~SubvolumeEvent()
        {
            ;
        }

        virtual void fire()
        {
            assert(proxy_ != NULL);
            sim_->reset_last_reactions();
            proxy_->fire(time_, coord_);
            update();
        }

        virtual void interrupt(Real const& t)
        {
            time_ = t;
            update();
        }

        void update()
        {
            const std::pair<Real, ReactionRuleProxyBase*>
                retval = sim_->draw_next_reaction(coord_);
            time_ += retval.first;
            proxy_ = retval.second;
        }

    protected:

        MesoscopicSimulator* sim_;
        coordinate_type coord_;

        ReactionRuleProxyBase* proxy_;
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
    virtual bool check_reaction() const
    {
        return last_reactions_.size() > 0;
    }

    std::vector<std::pair<ReactionRule, reaction_info_type> > last_reactions() const
    {
        return last_reactions_;
    }

    void add_last_reaction(const ReactionRule& rr, const reaction_info_type& ri)
    {
        last_reactions_.push_back(std::make_pair(rr, ri));
    }

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

    DiffusionProxy* create_diffusion_proxy(const Species& sp);

    void interrupt_all(const Real& t);
    std::pair<Real, ReactionRuleProxyBase*>
        draw_next_reaction(const coordinate_type& c);

    void increment_molecules(const Species& sp, const coordinate_type& c);
    void decrement_molecules(const Species& sp, const coordinate_type& c);
    void increment(const boost::shared_ptr<MesoscopicWorld::PoolBase>& pool, const coordinate_type& c);
    void decrement(const boost::shared_ptr<MesoscopicWorld::PoolBase>& pool, const coordinate_type& c);

protected:

    std::vector<std::pair<ReactionRule, reaction_info_type> > last_reactions_;

    boost::ptr_vector<ReactionRuleProxyBase> proxies_;
    boost::ptr_vector<ReactionRuleProxyBase>::size_type diffusion_proxy_offset_;

    EventScheduler scheduler_;
    std::vector<EventScheduler::identifier_type> event_ids_;
    coordinate_type interrupted_;
};

} // meso

} // ecell4

#endif /* ECELL4_MESO_MESOSCOPIC_SIMULATOR_HPP */
