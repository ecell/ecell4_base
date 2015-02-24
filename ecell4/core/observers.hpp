#ifndef __ECELL4_OBSERVER_HPP
#define __ECELL4_OBSERVER_HPP

#include "types.hpp"
#include "Space.hpp"
#include "Simulator.hpp"

#include <fstream>
#include <boost/format.hpp>


namespace ecell4
{

class Observer
{
public:

    Observer(const bool e)
        : every_(e)
    {
        ;
    }

    virtual ~Observer()
    {
        ; // do nothing
    }

    virtual const Real next_time() const
    {
        return inf;
    }

    virtual void initialize(const Space* space)
    {
        ;
    }

    virtual void finalize(const Space* space)
    {
        ;
    }

    virtual void reset()
    {
        ;
    }

    virtual void fire(const Simulator* sim, const Space* space) = 0;

    bool every()
    {
        return every_;
    }

private:

    const bool every_;
};

class FixedIntervalObserver
    : public Observer
{
public:

    typedef Observer base_type;

public:

    FixedIntervalObserver(const Real& dt)
        : base_type(false), t0_(0.0), dt_(dt), num_steps_(0), count_(0)
    {
        ;
    }

    virtual ~FixedIntervalObserver()
    {
        ;
    }

    const Real next_time() const
    {
        return t0_ + dt_ * count_;
    }

    const Integer num_steps() const
    {
        return num_steps_;
    }

    virtual void initialize(const Space* space)
    {
        t0_ = space->t();
        count_ = 0;
        // num_steps_ = 0;
    }

    virtual void fire(const Simulator* sim, const Space* space)
    {
        // tnext_ += dt_;
        ++count_;
        ++num_steps_;
    }

    virtual void reset()
    {
        num_steps_ = 0;
    }

protected:

    Real t0_, dt_;
    Integer num_steps_;
    Integer count_;
};

struct NumberLogger
{
    typedef std::vector<std::vector<Real> > data_container_type;
    typedef std::vector<Species> species_container_type;

    NumberLogger(const std::vector<std::string>& species)
    {
        targets.reserve(species.size());
        for (std::vector<std::string>::const_iterator i(species.begin());
            i != species.end(); ++i)
        {
            targets.push_back(Species(*i));
        }
    }

    ~NumberLogger()
    {
        ;
    }

    void initialize()
    {
        ;
    }

    void reset()
    {
        data.clear();
    }

    void log(const Space* space)
    {
        data_container_type::value_type tmp;
        tmp.push_back(space->t());
        for (species_container_type::const_iterator i(targets.begin());
            i != targets.end(); ++i)
        {
            tmp.push_back(space->get_value(*i));
            // tmp.push_back(space->num_molecules(*i));
        }
        data.push_back(tmp);
    }

    data_container_type data;
    species_container_type targets;
};

class FixedIntervalNumberObserver
    : public FixedIntervalObserver
{
public:

    typedef FixedIntervalObserver base_type;

public:

    FixedIntervalNumberObserver(const Real& dt, const std::vector<std::string>& species)
        : base_type(dt), logger_(species)
    {
        ;
    }

    virtual ~FixedIntervalNumberObserver()
    {
        ;
    }

    virtual void initialize(const Space* space)
    {
        base_type::initialize(space);
        logger_.initialize();
    }

    virtual void fire(const Simulator* sim, const Space* space)
    {
        logger_.log(space);
        base_type::fire(sim, space);
    }

    virtual void reset()
    {
        logger_.reset();
        base_type::reset();
    }

    NumberLogger::data_container_type data() const
    {
        return logger_.data;
    }

    NumberLogger::species_container_type targets() const
    {
        return logger_.targets;
    }

protected:

    NumberLogger logger_;
};

class NumberObserver
    : public Observer
{
public:

    typedef Observer base_type;

public:

    NumberObserver(const std::vector<std::string>& species)
        : base_type(true), logger_(species)
    {
        ;
    }

    virtual ~NumberObserver()
    {
        ;
    }

    virtual void initialize(const Space* space)
    {
        base_type::initialize(space);
        logger_.initialize();
        logger_.log(space);
    }

    virtual void finalize(const Space* space)
    {
        logger_.log(space);
        base_type::finalize(space);
    }

    virtual void fire(const Simulator* sim, const Space* space)
    {
        if (sim->last_reactions().size() > 0)
        {
            logger_.log(space);
        }
    }

    virtual void reset()
    {
        logger_.reset();
        base_type::reset();
    }

    NumberLogger::data_container_type data() const
    {
        return logger_.data;
    }

    NumberLogger::species_container_type targets() const
    {
        return logger_.targets;
    }

protected:

    NumberLogger logger_;
};

class TimingObserver
    : public Observer
{
public:

    typedef Observer base_type;

public:

    TimingObserver(const std::vector<Real>& t)
        : base_type(false), t_(t), num_steps_(0), count_(0)
    {
        ;
    }

    virtual ~TimingObserver()
    {
        ;
    }

    const Real next_time() const
    {
        if (count_ >= t_.size())
        {
            return inf;
        }
        return t_[count_];
    }

    const Integer num_steps() const
    {
        return num_steps_;
    }

    virtual void initialize(const Space* space)
    {
        ;
    }

    virtual void fire(const Simulator* sim, const Space* space)
    {
        ++num_steps_;
        ++count_;
    }

    virtual void reset()
    {
        num_steps_ = 0;
        count_ = 0;
    }

protected:

    std::vector<Real> t_;
    Integer num_steps_;
    Integer count_;
};

class TimingNumberObserver
    : public TimingObserver
{
public:

    typedef TimingObserver base_type;

public:

    TimingNumberObserver(
        const std::vector<Real>& t, const std::vector<std::string>& species)
        : base_type(t), logger_(species)
    {
        ;
    }

    virtual ~TimingNumberObserver()
    {
        ;
    }

    virtual void initialize(const Space* space)
    {
        base_type::initialize(space);
        logger_.initialize();
    }

    virtual void fire(const Simulator* sim, const Space* space)
    {
        logger_.log(space);
        base_type::fire(sim, space);
    }

    virtual void reset()
    {
        logger_.reset();
        base_type::reset();
    }

    NumberLogger::data_container_type data() const
    {
        return logger_.data;
    }

    NumberLogger::species_container_type targets() const
    {
        return logger_.targets;
    }

protected:

    NumberLogger logger_;
};

class FixedIntervalHDF5Observer
    : public FixedIntervalObserver
{
public:

    typedef FixedIntervalObserver base_type;

public:

    FixedIntervalHDF5Observer(const Real& dt, const std::string& filename)
        : base_type(dt), prefix_(filename)
    {
        ;
    }

    virtual ~FixedIntervalHDF5Observer()
    {
        ;
    }

    virtual void initialize(const Space* space)
    {
        base_type::initialize(space);
    }

    virtual void fire(const Simulator* sim, const Space* space)
    {
        space->save(filename());

        base_type::fire(sim, space);
    }

    const std::string filename() const
    {
        boost::format fmt(prefix_);

        if (fmt.expected_args() == 0)
        {
            return fmt.str();
        }
        else
        {
            return (fmt % num_steps()).str();
        }
    }

protected:

    std::string prefix_;
};

class FixedIntervalCSVObserver
    : public FixedIntervalObserver
{
public:

    typedef FixedIntervalObserver base_type;

public:

    FixedIntervalCSVObserver(const Real& dt, const std::string& filename)
        : base_type(dt), prefix_(filename)
    {
        ;
    }

    virtual ~FixedIntervalCSVObserver()
    {
        ;
    }

    virtual void initialize(const Space* space)
    {
        base_type::initialize(space);
    }

    virtual void fire(const Simulator* sim, const Space* space)
    {
        log(space);
        base_type::fire(sim, space);
    }

    void log(const Space* space)
    {
        typedef std::vector<std::pair<ParticleID, Particle> >
            particle_container_type;
        typedef utils::get_mapper_mf<Species::serial_type, unsigned int>::type
            serial_map_type;

        const particle_container_type particles(space->list_particles());
        serial_map_type serials;
        unsigned int cnt(0);

        std::ofstream ofs(filename().c_str(), std::ios::out);
        ofs << std::setprecision(17);
        ofs << "x,y,z,r,sid" << std::endl;
        for(particle_container_type::const_iterator i(particles.begin());
            i != particles.end(); ++i)
        {
            const Real3 pos((*i).second.position());
            const Real radius((*i).second.radius());

            unsigned int idx;
            serial_map_type::iterator
                j(serials.find((*i).second.species_serial()));
            if (j == serials.end())
            {
                idx = cnt;
                serials.insert(std::make_pair((*i).second.species_serial(), idx));
                ++cnt;
            }
            else
            {
                idx = (*j).second;
            }

            ofs << pos[0] << "," << pos[1] << "," << pos[2] << "," << radius
                << "," << idx << std::endl;
        }

        ofs.close();
    }

    const std::string filename() const
    {
        boost::format fmt(prefix_);

        if (fmt.expected_args() == 0)
        {
            return fmt.str();
        }
        else
        {
            return (fmt % num_steps()).str();
        }
    }

protected:

    std::string prefix_;
};

class FixedIntervalTrajectoryObserver
    : public FixedIntervalObserver
{
public:

    typedef FixedIntervalObserver base_type;

public:

    FixedIntervalTrajectoryObserver(
        const Real& dt, const std::vector<ParticleID>& pids,
        const bool& resolve_boundary = true)
        : base_type(dt), pids_(pids), resolve_boundary_(resolve_boundary)
    {
        ;
    }

    virtual ~FixedIntervalTrajectoryObserver()
    {
        ;
    }

    virtual void initialize(const Space* space)
    {
        base_type::initialize(space);
    }

    virtual void fire(const Simulator* sim, const Space* space)
    {
        base_type::fire(sim, space);

        const Real3 edge_lengths(space->edge_lengths());
        std::vector<std::vector<Real3> >::iterator j(trajectories_.begin());
        std::vector<Real3>::iterator k(strides_.begin());
        for (std::vector<ParticleID>::const_iterator i(pids_.begin());
            i != pids_.end(); ++i)
        {
            if (space->has_particle(*i))
            {
                Real3& stride(*k);
                Real3 pos(stride
                    + space->get_particle(*i).second.position());
                if (resolve_boundary_ && (*j).size() > 0)
                {
                    const Real3 prev((*j)[(*j).size() - 1]);
                    for (unsigned int dim(0); dim != 3; ++dim)
                    {
                        const Real L(edge_lengths[dim]);
                        if (pos[dim] - prev[dim] >= L * 0.5)
                        {
                            stride[dim] -= L;
                            pos[dim] -= L;
                        }
                        else if (pos[dim] - prev[dim] <= L * -0.5)
                        {
                            stride[dim] += L;
                            pos[dim] += L;
                        }
                    }
                }
                (*j).push_back(pos);
            }
            ++j;
            ++k;
        }
    }

    virtual void reset()
    {
        base_type::reset();
        trajectories_.clear();
        trajectories_.resize(pids_.size(), std::vector<Real3>());
        strides_.clear();
        strides_.resize(pids_.size(), Real3(0, 0, 0));
    }

    const std::vector<std::vector<Real3> >& data() const
    {
        return trajectories_;
    }

protected:

    std::vector<ParticleID> pids_;
    std::vector<std::vector<Real3> > trajectories_;
    std::vector<Real3> strides_;
    bool resolve_boundary_;
};

// class BioImagingObserver
//     : public Observer
// {
// public:
// 
//     typedef Observer base_type;
//     typedef utils::get_mapper_mf<Species::serial_type, unsigned int>::type
//         serial_map_type;
// 
// public:
// 
//     BioImagingObserver(const Real& dt, const Real& exposure_time, const Integer& num_div, const Real& voxel_radius, const Real& scale)
//         : base_type(false), t0_(0.0), dt_(dt), exposure_time_(exposure_time),
//         num_div_(num_div), voxel_radius_(voxel_radius), scale_(scale), num_steps_(0)
//     {
//         ;
//     }
// 
//     virtual ~BioImagingObserver()
//     {
//         ;
//     }
// 
//     const Real next_time() const
//     {
//         if (num_div_ > 1)
//         {
//             const Real offset(t0_ + dt_ * static_cast<Real>(num_steps_ / num_div_));
//             return offset + (exposure_time_ / num_div_) * (num_steps_ % num_div_);
//         }
//         else
//         {
//             return t0_ + dt_ * num_steps_;
//         }
//     }
// 
//     const Integer num_steps() const
//     {
//         return num_steps_;
//     }
// 
//     virtual void initialize(const Space* space)
//     {
//         t0_ = space->t();
//         serials_.clear();
//         num_steps_ = 0;
//     }
// 
//     virtual void reset()
//     {
//         ;
//     }
// 
//     virtual void fire(const Simulator* sim, const Space* space)
//     {
//         log(space);
//         ++num_steps_;
//     }
// 
//     void log(const Space* space)
//     {
//         typedef std::vector<std::pair<ParticleID, Particle> >
//             particle_container_type;
//         const particle_container_type particles(space->list_particles());
// 
//         unsigned int cnt(serials_.size());
// 
//         const Real t(space->t());
// 
//         std::ofstream ofs(filename().c_str(), std::ios::out);
//         ofs << std::setprecision(17);
//         for(particle_container_type::const_iterator i(particles.begin());
//             i != particles.end(); ++i)
//         {
//             const Real3 pos((*i).second.position());
//             const Real radius((*i).second.radius());
//             const ParticleID& pid((*i).first);
// 
//             unsigned int idx;
//             serial_map_type::iterator
//                 j(serials_.find((*i).second.species_serial()));
//             if (j == serials_.end())
//             {
//                 idx = cnt;
//                 serials_.insert(std::make_pair((*i).second.species_serial(), idx));
//                 ++cnt;
//             }
//             else
//             {
//                 idx = (*j).second;
//             }
// 
//             ofs << t
//                 << "," << pos[0] * scale_
//                 << "," << pos[1] * scale_
//                 << "," << pos[2] * scale_
//                 << "," << radius * scale_
//                 << ",\"(" << pid.serial() << ",0)\""
//                 << ",\"(" << idx << ",0)\"" << std::endl;
//         }
// 
//         ofs.close();
// 
//         write_header(space);
//     }
// 
//     void write_header(const Space* space)
//     {
//         std::ofstream ofs("pt-input.csv", std::ios::out);
// 
//         if (num_div_ > 1)
//         {
//             ofs << exposure_time_ / num_div_;
//         }
//         else
//         {
//             ofs << exposure_time_;
//         }
// 
//         const Real3 edge_lengths(space->edge_lengths());
//         ofs << "," << edge_lengths[1] / (voxel_radius_ * 2)
//             << "," << edge_lengths[2] / (voxel_radius_ * 2)
//             << "," << edge_lengths[0] / (voxel_radius_ * 2)
//             << "," << voxel_radius_ * scale_;
// 
//         std::vector<Species::serial_type> species(serials_.size());
//         for (serial_map_type::const_iterator i(serials_.begin());
//             i != serials_.end(); ++i)
//         {
//             species[(*i).second] = (*i).first;
//         }
//         for (std::vector<Species::serial_type>::const_iterator i(species.begin());
//             i != species.end(); ++i)
//         {
//             ofs << ",[/:" << (*i) << "]=" << voxel_radius_ * scale_;
//         }
//         ofs << std::endl;
//         ofs.close();
//     }
// 
//     const std::string filename() const
//     {
//         const Integer i(num_steps_ / num_div_);
//         const Integer j(num_steps_ % num_div_);
// 
//         boost::format fmt("pt-%09d.%03d.csv");
//         const std::string fname((fmt % i % j).str());
//         return fname;
//     }
// 
// protected:
// 
//     Real t0_, dt_, exposure_time_;
//     Integer num_div_;
//     Real voxel_radius_, scale_;
//     Integer num_steps_;
// 
//     serial_map_type serials_;
// };

} // ecell4

#endif /* __ECELL4_OBSEVER_HPP */
