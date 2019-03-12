#ifndef ECELL4_SGFRD_BD_SIMULATOR
#define ECELL4_SGFRD_BD_SIMULATOR

#include <ecell4/core/SimulatorBase.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/geometry.hpp>
#include <ecell4/sgfrd/ReactionInfo.hpp>
#include <ecell4/sgfrd/SGFRDWorld.hpp>
#include <ecell4/sgfrd/BDPropagator.hpp>
#include <ecell4/sgfrd/tracer.hpp>

namespace ecell4
{
namespace sgfrd
{

class BDSimulator :
    public ecell4::SimulatorBase<SGFRDWorld, ecell4::Model>
{
  public:
    typedef BDSimulator self_type;

    // polygon
    typedef ecell4::Polygon  polygon_type;
    typedef polygon_type::FaceID   FaceID;
    typedef polygon_type::EdgeID   EdgeID;
    typedef polygon_type::VertexID VertexID;

    // world & particles
    typedef ecell4::SimulatorBase<SGFRDWorld, ecell4::Model> base_type;
    typedef base_type::world_type world_type;
    typedef base_type::model_type model_type;
    typedef std::pair<ParticleID, Particle> particle_id_pair_type;
    typedef boost::tuple<ParticleID, Particle, FaceID> pid_p_fid_tuple_type;

    // reaction
    typedef ecell4::ReactionRule           reaction_rule_type;
    typedef ReactionInfo                   reaction_info_type;
    typedef std::pair<reaction_rule_type, reaction_info_type> reaction_log_type;
    typedef std::vector<reaction_log_type> reaction_archive_type;

    BDSimulator(const boost::shared_ptr<world_type>& world,
                const boost::shared_ptr<model_type>& model,
                Real bd_dt_factor = 1e-5, Real reaction_length = 1e-3,
                const std::string& trace_fname = "bd_trace.log")
        : base_type(world, model), dt_(bd_dt_factor/*TODO*/),
          bd_dt_factor_(bd_dt_factor), reaction_length_(reaction_length),
          rng_(*(world->rng())), tracer_(trace_fname)
    {}

    BDSimulator(boost::shared_ptr<world_type> world,
                Real bd_dt_factor = 1e-5, Real reaction_length = 1e-3,
                const std::string& trace_fname = "bd_trace.log")
        : base_type(world), dt_(bd_dt_factor/*TODO*/),
          bd_dt_factor_(bd_dt_factor), reaction_length_(reaction_length),
          rng_(*(world->rng())), tracer_(trace_fname)
    {}

    void initialize()
    {
        return;
    }
    void finalize()
    {
        return;
    }
    void finalize(const Real t)
    {
        SGFRD_SCOPE(us, finalize, tracer_);

        this->world_->set_t(t);
        BDPropagator<world_type, volume_clearer> propagator(
                *(this->model_), *(this->world_),
                *(this->world_->polygon()), this->rng_,
                t - this->world_->t(), this->reaction_length_,
                this->last_reactions_,
                volume_clearer(*(this->world_), this->tracer_));

        while(propagator())
        {
            // do nothing
        }
        return;
    }

    void step()
    {
        SGFRD_SCOPE(us, step, tracer_);

        const Real next_time = this->world_->t() + this->dt_;
        this->world_->set_t(next_time);
        SGFRD_TRACE(tracer_.write("now t = %1%", this->world_->t()))
        SGFRD_TRACE(tracer_.write("dt = %1%", this->dt_))

        BDPropagator<world_type, volume_clearer> propagator(
                *(this->model_), *(this->world_),
                *(this->world_->polygon()), this->rng_,
                this->dt_, this->reaction_length_,
                this->last_reactions_,
                volume_clearer(*(this->world_), this->tracer_));

        SGFRD_TRACE(tracer_.write("propagating..."))
        while(propagator())
        {
            // do nothing
        }
        SGFRD_TRACE(tracer_.write("...done!"))
        return;
    }
    bool step(const Real& upto)
    {
        this->step();
        return this->world_->t() < upto;
    }

    Real dt() const {return dt_;}
    Real reaction_length() const {return reaction_length_;}

    bool check_reaction() const {return last_reactions_.size() > 0;}

    std::vector<std::pair<ReactionRule, reaction_info_type> > const&
    last_reactions() const {return last_reactions_;}

    Real next_event_time() const
    {
        return this->world_->t() + this->dt_;
    }

    struct volume_clearer
    {
        volume_clearer(world_type const& wld, tracer& tr)
            : world_(wld), tracer_(tr)
        {}

        bool operator()(const Particle& p, const FaceID& fid) const
        {
            return world_.check_no_overlap(std::make_pair(p.position(), fid),
                                           p.radius());
        }
        bool operator()(const Particle& p, const FaceID& fid,
                        const ParticleID& ignore) const
        {
            return world_.check_no_overlap(std::make_pair(p.position(), fid),
                                           p.radius(), ignore);
        }
        bool operator()(const Particle& p, const FaceID& fid,
                        const ParticleID& ign1, const ParticleID& ign2) const
        {
            return world_.check_no_overlap(std::make_pair(p.position(), fid),
                                           p.radius(), ign1, ign2);
        }

        tracer& access_tracer() {return this->tracer_;}

        const world_type& world_;
        tracer& tracer_;
    };

  private:
    // from SimulatorBase
    // boost::shared_ptr<model_type> model_;
    // boost::shared_ptr<world_type> world_;
    // Integer num_steps_;

    Real dt_;
    Real bd_dt_factor_;
    Real reaction_length_;
    ecell4::RandomNumberGenerator& rng_;
    std::vector<std::pair<reaction_rule_type, reaction_info_type> > last_reactions_;
    mutable tracer tracer_;
};

} // sgfrd
} // ecell4
#endif// ECELL4_SGFRD_BD_SIMULATOR
