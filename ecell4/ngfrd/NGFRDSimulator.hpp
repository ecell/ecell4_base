#ifndef ECELL4_NGFRD_NGFRDSIMULATOR_HPP
#define ECELL4_NGFRD_NGFRDSIMULATOR_HPP

#include <boost/format.hpp>
#include <boost/none_t.hpp>
#include <boost/optional.hpp>
#include <boost/variant.hpp>
#include <boost/container/small_vector.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>

#include <ecell4/ngfrd/ShellID.hpp>
#include <ecell4/ngfrd/Shell.hpp>
#include <ecell4/ngfrd/DomainID.hpp>
#include <ecell4/ngfrd/Domain.hpp>
#include <ecell4/ngfrd/NGFRDEvent.hpp>
#include <ecell4/ngfrd/NGFRDWorld.hpp>

#include <greens_functions/PairGreensFunction.hpp>
#include <greens_functions/GreensFunction3DRadAbs.hpp>
#include <greens_functions/GreensFunction3DRadInf.hpp>
#include <greens_functions/GreensFunction3DAbsSym.hpp>
#include <greens_functions/GreensFunction3DAbs.hpp>
#include <greens_functions/GreensFunction3D.hpp>

namespace ecell4
{
namespace ngfrd
{

class NGFRDSimulator final: public SimulatorBase<NGFRDWorld>
{
public:
    static constexpr Real SAFETY              = 1.0 + 1e-5;
    static constexpr Real SINGLE_SHELL_FACTOR = 0.1;
    static constexpr Real DEFAULT_DT_FACTOR   = 1e-5;
    static constexpr Real CUTOFF_FACTOR       = 5.6;

public:

    using base_type       = SimulatorBase<NGFRDWorld>;
    using model_type      = typename base_type::model_type;
    using world_type      = typename base_type::world_type;
    using event_type      = NGFRDEvent;
    using scheduler_type  = EventSchedulerBase<event_type>;
    using event_id_type   = typename scheduler_type::identifier_type;

    NGFRDSimulator(const std::shared_ptr<world_type>& world,
                   const std::shared_ptr<model_type>& model,
                   const Real bd_dt_factor_3D = 1e-5,
                   const Real bd_dt_factor_2D = 1e-3,
                   const Real reaction_length = 1e-1,
                   const std::size_t max_retry_moves = 1)
        : base_type(world, model), dt_factor_3D_(bd_dt_factor_3D),
          dt_factor_2D_(bd_dt_factor_2D), reaction_length_(reaction_length),
          max_retry_(max_retry_moves),
          shells_(world->edge_lengths(), world->polygon_ptr()),
          is_uninitialized_(true)
    {}
    ~NGFRDSimulator() override = default;

    void initialize() override
    {
        // --------------------------------------------------------------------
        // clear everything

        scheduler_.clear();
        domains_  .clear();
        shells_   .clear();
        shells_.reset_boundary(this->world_.edge_lengths());

        // --------------------------------------------------------------------
        // form domain for all particles

        for(const auto& pp : world_->particles())
        {
            form_domain(pp.first, pp.second);
        }

        // --------------------------------------------------------------------
        // form birth domain if needed

        // TODO!

        return;
    }
    void step() override
    {
        if(this->is_uninitialized_)
        {
            this->initialize();
        }
        this->step_unchecked();
        return;
    }

    bool step(const Real& upto) override
    {
        if(this->is_uninitialized_)
        {
            this->initialize();
        }

        if(upto <= this->t())
        {
            return false;
        }

        if(scheduler_.next_time() <= upto)
        {
            this->step_unchecked();
            return true;
        }
        else
        {
            this->set_t(upto);
            this->finalize();
            return false;
        }
    }
    void finalize()
    {
        std::vector<domain_id_type> non_singles;

        // Single does not stick out from shell when it is bursted.
        for(const auto& eidp: scheduler_.events())
        {
            const auto& eid = eidp.first;
            const auto& ev  = eidp.second;
            const auto& dom = this->get_domain(ev.domain_id());

            if(dom.multiplicity() != 1)
            {
                non_singles.push_back(ev.domain_id());
                continue;
            }
            this->burst_domain(dom);
        }

        // then burst non-single domains
        for(const auto& did : non_singles)
        {
            this->burst_domain(this->get_domain(did));
        }
        this->dt_ = 0.0;
        this->is_uninitialized_ = true;

        scheduler_.clear();
        return;
    }

    Real next_time() const override
    {
        return scheduler_.next_time();
    }
    Real dt() const override
    {
        return scheduler_.next_time() - world_->t();
    }

    std::shared_ptr<RandomNumberGenerator> const& rng() const noexcept
    {
        return this->world_->rng();
    }

private:

    Domain const& get_domain(const DomainID& did) const {return domains_.at(did).second;}
    Domain&       get_domain(const DomainID& did)       {return domains_.at(did).second;}

    void step_unchecked()
    {
        this->num_steps_ += 1;

        if(scheduler_.size() == 0)
        {
            this->set_t(scheduler_.next_time());
            return;
        }

        const auto eidp = scheduler_.pop();
        this->set_t(eidp.second->time());

        const auto fired = fire_event(eidp.second);
        for(const auto& pidp : fired)
        {
            this->form_domain(pidp.first, pidp.second);
        }

        const auto next_time = scheduler_.top().second->time();
        this->dt_ = next_time - this->t();

        if(this->dt_ == 0.0)
        {
            this->zero_step_count_ += 1;
        }
        return;
    }

    void form_domain(const ParticleID& pid, const Particle& p)
    {
        if(const auto fid = this->world_.on_which_face(pid))
        {
            this->form_domain_2D(pid, p, *fid);
        }
        else
        {
            this->form_domain_3D(pid, p);
        }
    }
    void form_domain_2D(const ParticleID& pid, const Particle& p);
    void form_domain_3D(const ParticleID& pid, const Particle& p);

    void fire_event(const NGFRDEvent& ev)
    {
        return fire_domain(ev.domain_id(), this->get_domain(ev.domain_id()));
    }

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
    fire_domain(const DomainID& did, const Domain& dom)
    {
        switch(dom.kind())
        {
            // TODO: add more
            case Domain::DomainKind::Multi:
            {
                return this->fire_multi(did, dom.as_multi());
            }
            default:
            {
                throw_exception<NotImplemented>("NGFRD::fire_domain: unknown "
                        "domain kind (", static_cast<int>(dom.kind()), ").");
            }
        }
    }

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
    fire_multi(const DomainID& did, const MultiDomain& dom);

private:

    // inherited from SimulatorBase
//     std::shared_ptr<world_type> world_;
//     std::shared_ptr<model_type> model_;
//     Integer num_steps_;

    Real dt_factor_3D_;
    Real dt_factor_2D_;
    Real reaction_length_;
    std::size_t max_retry_;

    SerialIDGenerator<DomainID> didgen_;
    std::unordered_map<DomainID, std::pair<event_id_type, Domain>> domains_;

    SerialIDGenerator<ShellId> sidgen_;
    ShellContainer             shells_;

    scheduler_type scheduler_;
    Real dt_;
    bool is_uninitialized_;

    // -------------------------------------------------------------------------
    // statistics

    std::size_t rejected_moves_;
    std::size_t zero_step_count_;
//     std::array<std::size_t, SingleDomain::EventKinds> single_step_count_;
//     std::array<std::size_t, PairDomain::EventKinds>   pair_step_count_;
//     std::array<std::size_t, MultiDomain::EventKinds>  multi_step_count_;
//     std::array<std::size_t, NUM_DOMAIN_KINDS>         domain_count_per_type_;
};
#undef CHECK

} // egfrd
} // ecell4
#endif /* EGFRDSIMULATOR_HPP */
