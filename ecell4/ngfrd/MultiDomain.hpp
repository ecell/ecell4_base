#ifndef ECELL4_NGFRD_MULTI_DOMAIN_HPP
#define ECELL4_NGFRD_MULTI_DOMAIN_HPP
#include <ecell4/ngfrd/ShellID.hpp>
#include <ecell4/ngfrd/BDPropagator.hpp>
#include <ecell4/ngfrd/NGFRDWorld.hpp>
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/Particle.hpp>
#include <ecell4/core/Model.hpp>

#include <boost/container/small_vector.hpp>

namespace ecell4
{
namespace ngfrd
{

class NGFRDSimulator; // forwad declaration
class NGFRDWorld;     // ditto

class MultiDomain
{
  public:

    enum class EventKind : std::uint8_t
    {
        None,
        Escape,
        Reaction,
    };
    using reaction_log_type = std::pair<ReactionRule, ReactionInfo>;
    using reaction_log_container_type =
        boost::container::small_vector<reaction_log_type, 2>;

    // To avoid an allocation and to gain memory locality, we use small_vector
    // of Identifiers. Essentially Multi can contain arbitrary number of shells
    // and particles, so we need to use small_vector instead of static_vector.
    // We need to find the sweet spot between the risk of memory loss caused by
    // too large statically allocated region in the small_vector and the benefit
    // from avoiding allocation and locality cost. Now we chose 8, but this can
    // be changed.
    //     The following makes this class relatively bigger (it contains static
    // array of length 8). It is inefficient to copy/move this class with
    // others. This means that the domain struct should not be included in the
    // Event class. Event should contain DomainID, and simulator shopuld map
    // DomainID to the actual Domain in a container.
    using shell_id_container_type =
        boost::container::small_vector<ShellID, 8>;
    using particle_id_container_type =
        boost::container::small_vector<ParticleID, 8>;

  public:

    MultiDomain(const Real begin_t,
                const Real dt_factor_3D = 1.0,
                const Real dt_factor_2D = 0.01,
                const Real rl_factor    = 0.1) noexcept
        : kind_(None), begin_time_(begin_t),
          dt_(-1.0), dt_factor_3D_(dt_factor_3D), dt_factor_2D_(dt_factor_2D)
          reaction_length_factor_(rl_factor), reaction_length_(-1.0)
    {}
    ~MultiDomain() = default;

    // if inserted (i.e. not already inserted), returns true.
    bool add_particle(const ParticleID& pid)
    {
        if(std::find(particle_ids_.begin(), particle_ids_.end(), pid) !=
                     particle_ids_.end()) // found?
        {
            return false;
        }
        this->particle_ids_.push_back(pid);
        return true;
    }
    bool add_shell(const ShellID& sid)
    {
        if(std::find(shell_ids_.begin(), shell_ids_.end(), sid) !=
                     shell_ids_.end()) // found?
        {
            return false;
        }
        this->shell_ids_.push_back(sid);
        return true;
    }

    void step(const Model& model, NGFRDSimulator& sim, NGFRDWorld& world)
    {
        this->step(model, sim, world, this->dt_);
        return ;
    }

    void step(const Model& model, NGFRDSimulator& sim, NGFRDWorld& world, const Real dt);

    EventKind& eventkind()       noexcept {return this->kind_;}
    EventKind  eventkind() const noexcept {return this->kind_;}

    Real& dt()                    noexcept {return this->dt_;}
    Real  dt()              const noexcept {return this->dt_;}
    Real& begin_time()            noexcept {return this->begin_time_;}
    Real  begin_time()      const noexcept {return this->begin_time_;}
    Real& reaction_length()       noexcept {return this->reaction_length_;}
    Real  reaction_length() const noexcept {return this->reaction_length_;}

    shell_id_container_type&          shell_ids()          noexcept {return shells_;}
    shell_id_container_type    const& shell_ids()    const noexcept {return shells_;}
    particle_id_container_type&       particle_ids()       noexcept {return particles_;}
    particle_id_container_type const& particle_ids() const noexcept {return particles_;}

    std::size_t num_shells()   const noexcept {return shells_.size();}
    std::size_t multiplicity() const noexcept {return particles_.size();}

    reaction_log_container_type const& last_reactions() const {return last_reactions_;}

    void determine_parameters(const Model& model, const NGFRDWorld& world)
    {
        this->determine_reaction_length(model, world);
        this->determine_delta_t(model, world);
        return;
    }

  private:

    void determine_reaction_length(const Model& model, const NGFRDWorld& world)
    {
        Real min_radius = std::numeric_limits<Real>::max();
        for(const auto& pid : this->particle_ids())
        {
            min_radius = std::min(world.get_particle(pid).second.radius(), min_radius);
        }
        this->reaction_length_ = this->reaction_length_factor_ * min_radius;
        return;
    }
    void determine_delta_t(const Model& model, const NGFRDWorld& world)
    {
        const Real dt_2D = this->determine_delta_t_2D(model, world);
        const Real dt_3D = this->determine_delta_t_3D(model, world);
        this->dt_ = std::min(dt_2D, dt_3D);
        return;
    }

    Real determine_delta_t_3D(const Model& model, const NGFRDWorld& world)
    {
        // check birth events
        Real birth_prob = 0;
        for(const auto& rule : model.reaction_rules())
        {
            if(!rule.reactants().empty())
            {
                continue;
            }
            birth_prob += rule.k();
        }
        birth_prob *= world.volume();

        std::vector<ParticleID> p3D;
        for(const auto& pid : this->particle_ids())
        {
            if(world.on_which_face(pid).has_value())
            {
                p3D.push_back(pid);
            }
        }

        if(p3D.empty() && birth_prob == 0.0)
        {
            return std::numeric_limits<Real>::infinity(); // anything is okay.
        }

        Real D_max(0.0), radius_min(std::numeric_limits<Real>::max());
        for(const auto& pid : p3D)
        {
            const auto& species = world.get_particle(pid).species();
            const auto  molinfo = world.get_molecule_info(species);
            D_max      = std::max(molinfo.D,      D_max);
            radius_min = std::min(molinfo.radius, radius_min);
        }
        const Real dt = radius_min * radius_min * 4 / (D_max * 2);

        if(birth_prob == 0.0)
        {
            return dt;
        }
        else
        {
            return std::min(dt, 1.0 / prob);
        }
    }

    Real determine_delta_t_2D(const Model& model, const NGFRDWorld& world)
    {
        // it assumes this->determine_reaction_length() has already been called
        assert(this->reaction_length_ > 0.0);

        // check birth events
        Real birth_prob = 0;
        for(const auto& rule : model.reaction_rules())
        {
            if(!rule.reactants().empty())
            {
                continue;
            }
            birth_prob += rule.k();
        }
        birth_prob *= world.volume();

        std::vector<ParticleID> p2D;
        for(const auto& pid : this->particle_ids())
        {
            if(world.on_which_face(pid).has_value())
            {
                p2D.push_back(pid);
            }
        }

        // there could be a birth reaction...
        if(p2D.empty() && birth_prob == 0.0)
        {
            return std::numeric_limits<Real>::infinity(); // anything is okay
        }

        // find the maximum D in 2D particles
        std::vector<Species> sps;
        sps.reserve(p2D.size());
        Real D_max = -std::numeric_limits<Real>::max();
        for(const auto& pid : p2D)
        {
            const auto& species = world.get_particle(pid).species();
            const auto  molinfo = world.get_molecule_info(species);
            D_max = std::max(molinfo.D, D_max);

            if(std::find(sps.begin(), sps.end(), species) != sps.end())
            {
                sps.push_back(species);
            }
        }

        // find the maximum reaction rate in 2 particle reaciton
        Real k_max = 0.0;
        for(std::size_t i=0; i<sps.size(); ++i)
        {
            const auto& sp1 = sps.at(i);
            // start from j = i. there can be a reaction like "A+A -> B".
            for(std::size_t j=i; j<sps.size(); ++j)
            {
                const auto& sp2  = sps.at(j);
                for(auto&& rule : model_.query_reaction_rules(sp1, sp2))
                {
                    k_max = std::max(k_max, rule.k());
                }
            }
        }

        // calculate dt
        const Real P_max = 0.05; // upper limit for reaction probability
        const Real delta = this->reaction_length_;

        const Real upper_limit_from_D = delta * delta / D_max;
        const Real upper_limit_from_k = (k_max > 0.0) ? (P_max * delta / k_max) :
                                        std::numeric_limits<Real>::infinity();

        return dt_factor_2D * std::min(upper_limit_from_D, upper_limit_from_k);
    }

  private:

    EventKind kind_;
    Real begin_time_;
    Real dt_;
    Real dt_factor_3D_;
    Real dt_factor_2D_;
    Real reaction_length_;
    Real reaction_length_factor_;
    shell_id_container_type     shell_ids_;
    particle_id_container_type  particle_ids_;
    reaction_log_container_type last_reactions_;
};

} // sgfrd
} // ecell4
#endif //ECELL4_NGFRD_MULTI_DOMAIN_HPP
