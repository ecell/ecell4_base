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
        assert(this->diagnosis());
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

    bool diagnosis()
    {
        SGFRD_SCOPE(us, diagnosis, tracer_);
        bool result = true;
        auto particles = this->world_->list_particles();

        ParticleID pid; Particle p;
        for(const auto& pidp : particles)
        {
            std::tie(pid, p) = pidp;
            const FaceID fid = this->world_->get_face_id(pid);
            std::pair<Real3, FaceID> pos = std::make_pair(p.position(), fid);

            ParticleID _pid; Particle _p;
            for(const auto& _pidp : particles)
            {
                std::tie(_pid, _p) = _pidp;
                if(pid == _pid) {continue;}
                const FaceID _fid = this->world_->get_face_id(_pid);
                const Real dist = this->world_->distance(pos,
                                      std::make_pair(_p.position(), _fid));
                if(dist < p.radius() + _p.radius())
                {
                    result = false;
                    std::cerr << "ERROR: particle " << pid << " and " << _pid
                              << "overlaps!\n";
                    std::cerr << "     : distance = " << dist << " < sum of radii = "
                              << p.radius() + _p.radius() << '\n';
                    std::cerr << "     : particle " << pid << " has radius "
                              << p.radius() << " at " << p.position() << " on "
                              << fid << '\n';
                    std::cerr << "     : particle " << _pid << " has radius "
                              << _p.radius() << " at " << _p.position() << " on "
                              << _fid << '\n';
                }
            }

            if(!this->world_->polygon()->is_inside_of_boundary(p.position()))
            {
                std::cerr << "ERROR: particle " << pid << " is outside of the boundary!\n";
                std::cerr << "     : position = " << p.position()
                          << ", boundary = " << this->world_->polygon()->edge_lengths() << "\n";
                result = false;
            }

            const Triangle&   tri  = this->world_->polygon()->triangle_at(fid);
            const Barycentric bary = ::ecell4::to_barycentric(p.position(), tri);
            if(!is_inside(bary))
            {
                std::cerr << "ERROR: particle " << pid << " is not on the face " << fid << "\n";
                std::cerr << "     : position    = " << p.position() << ", face = " << tri << "\n";
                std::cerr << "     : barycentric = " << bary << ", face = " << tri << "\n";
                result = false;
            }
        }
        return result;
    }

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
