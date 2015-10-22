#include "observers.hpp"


namespace ecell4
{

const Real Observer::next_time() const
{
    return inf;
}

void Observer::initialize(const Space* space)
{
    ;
}

void Observer::finalize(const Space* space)
{
    ;
}

void Observer::reset()
{
    ;
}

const Real FixedIntervalObserver::next_time() const
{
    return t0_ + dt_ * count_;
}

const Integer FixedIntervalObserver::num_steps() const
{
    return num_steps_;
}

const Integer FixedIntervalObserver::count() const
{
    return count_;
}

void FixedIntervalObserver::initialize(const Space* space)
{
    if (count_ == 0)
    {
        t0_ = space->t();
    }
    else
    {
        while (next_time() < space->t())
        {
            ++count_;
        }
    }
}

bool FixedIntervalObserver::fire(const Simulator* sim, const Space* space)
{
    ++num_steps_;
    ++count_;
    return true;
}

void FixedIntervalObserver::reset()
{
    num_steps_ = 0;
    count_ = 0;
    t0_ = 0.0; //DUMMY
}

void FixedIntervalNumberObserver::initialize(const Space* space)
{
    base_type::initialize(space);
    logger_.initialize();
}

bool FixedIntervalNumberObserver::fire(const Simulator* sim, const Space* space)
{
    logger_.log(space);
    return base_type::fire(sim, space);
}

void FixedIntervalNumberObserver::reset()
{
    logger_.reset();
    base_type::reset();
}

NumberLogger::data_container_type FixedIntervalNumberObserver::data() const
{
    return logger_.data;
}

NumberLogger::species_container_type FixedIntervalNumberObserver::targets() const
{
    return logger_.targets;
}

void NumberObserver::initialize(const Space* space)
{
    base_type::initialize(space);
    logger_.initialize();
    logger_.log(space);
}

void NumberObserver::finalize(const Space* space)
{
    if (logger_.data.size() == 0 || logger_.data.back()[0] != space->t())
    {
        logger_.log(space);
    }
    base_type::finalize(space);
}

bool NumberObserver::fire(const Simulator* sim, const Space* space)
{
    if (sim->check_reaction())
    {
        logger_.log(space);
        ++num_steps_;
    }
    return true;
}

void NumberObserver::reset()
{
    num_steps_ = 0;
    logger_.reset();
    base_type::reset();
}

const Integer NumberObserver::num_steps() const
{
    return num_steps_;
}

NumberLogger::data_container_type NumberObserver::data() const
{
    return logger_.data;
}

NumberLogger::species_container_type NumberObserver::targets() const
{
    return logger_.targets;
}

const Real TimingObserver::next_time() const
{
    if (count_ >= static_cast<Integer>(t_.size()))
    {
        return inf;
    }
    return t_[count_];
}

void TimingObserver::initialize(const Space* space)
{
    while (next_time() < space->t())
    {
        ++count_;
    }
}

bool TimingObserver::fire(const Simulator* sim, const Space* space)
{
    ++num_steps_;
    ++count_;
    return true;
}

void TimingObserver::reset()
{
    num_steps_ = 0;
    count_ = 0;
}

void TimingNumberObserver::initialize(const Space* space)
{
    base_type::initialize(space);
    logger_.initialize();
}

bool TimingNumberObserver::fire(const Simulator* sim, const Space* space)
{
    logger_.log(space);
    return base_type::fire(sim, space);
}

void TimingNumberObserver::reset()
{
    logger_.reset();
    base_type::reset();
}

NumberLogger::data_container_type TimingNumberObserver::data() const
{
    return logger_.data;
}

NumberLogger::species_container_type TimingNumberObserver::targets() const
{
    return logger_.targets;
}

void FixedIntervalHDF5Observer::initialize(const Space* space)
{
    base_type::initialize(space);
}

bool FixedIntervalHDF5Observer::fire(const Simulator* sim, const Space* space)
{
    space->save(filename());

    return base_type::fire(sim, space);
}

const std::string FixedIntervalHDF5Observer::filename() const
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

void FixedIntervalCSVObserver::initialize(const Space* space)
{
    base_type::initialize(space);
}

bool FixedIntervalCSVObserver::fire(const Simulator* sim, const Space* space)
{
    log(space);
    return base_type::fire(sim, space);
}

void FixedIntervalCSVObserver::write_particles(
    std::ofstream& ofs, const particle_container_type& particles,
    const Species::serial_type label)
{
    for(particle_container_type::const_iterator i(particles.begin());
        i != particles.end(); ++i)
    {
        const Real3 pos((*i).second.position());
        const Real radius((*i).second.radius());
        const Species::serial_type serial(
            label == "" ? (*i).second.species_serial() : label);

        unsigned int idx;
        serial_map_type::iterator j(serials_.find(serial));
        if (j == serials_.end())
        {
            idx = serials_.size();
            serials_.insert(std::make_pair(serial, idx));
        }
        else
        {
            idx = (*j).second;
        }

        ofs << pos[0] << "," << pos[1] << "," << pos[2] << "," << radius
            << "," << idx << std::endl;
    }
}

void FixedIntervalCSVObserver::log(const Space* space)
{
    std::ofstream ofs(filename().c_str(), std::ios::out);
    ofs << std::setprecision(17);
    ofs << "x,y,z,r,sid" << std::endl;

    if (species_.size() == 0)
    {
        const particle_container_type particles(space->list_particles());
        write_particles(ofs, particles);
    }
    else
    {
        for (std::vector<std::string>::const_iterator i(species_.begin());
            i != species_.end(); ++i)
        {
            const Species sp(*i);
            const particle_container_type particles(space->list_particles(sp));
            write_particles(ofs, particles, *i);
        }
    }

    ofs.close();
}

const std::string FixedIntervalCSVObserver::filename() const
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

void FixedIntervalCSVObserver::reset()
{
    serials_.clear();
    base_type::reset();
}

void FixedIntervalTrajectoryObserver::initialize(const Space* space)
{
    base_type::initialize(space);

    typedef std::vector<std::pair<ParticleID, Particle> > particle_id_pairs;
    if (pids_.size() == 0)
    {
        particle_id_pairs const particles(space->list_particles());
        pids_.reserve(particles.size());
        trajectories_.resize(particles.size());
        strides_.resize(particles.size());
        for (particle_id_pairs::const_iterator i(particles.begin());
            i != particles.end(); ++i)
        {
            pids_.push_back((*i).first);
        }
        assert(pids_.size() == particles.size());
    }
}

bool FixedIntervalTrajectoryObserver::fire(const Simulator* sim, const Space* space)
{
    const bool retval = base_type::fire(sim, space);
    t_.push_back(space->t());

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

    return retval;
}

void FixedIntervalTrajectoryObserver::reset()
{
    base_type::reset();

    trajectories_.clear();
    trajectories_.resize(pids_.size(), std::vector<Real3>());
    strides_.clear();
    strides_.resize(pids_.size(), Real3(0, 0, 0));
    t_.clear();
}

const std::vector<std::vector<Real3> >& FixedIntervalTrajectoryObserver::data() const
{
    return trajectories_;
}

const Integer FixedIntervalTrajectoryObserver::num_tracers() const
{
    return pids_.size();
}

const std::vector<Real>& FixedIntervalTrajectoryObserver::t() const
{
    return t_;
}

void TimeoutObserver::initialize(const Space* space)
{
    base_type::initialize(space);
    duration_ = 0.0;
    time(&tstart_);
}

void TimeoutObserver::finalize(const Space* space)
{
    base_type::finalize(space);
    acc_ += duration_;
}

bool TimeoutObserver::fire(const Simulator* sim, const Space* space)
{
    time_t tnow;
    time(&tnow);
    duration_ = difftime(tnow, tstart_);
    if (duration_ >= interval_)
    {
        return false;
    }
    return true;
}

void TimeoutObserver::reset()
{
    base_type::reset();
    duration_ = 0.0;
    acc_ = 0.0;
    time(&tstart_);
}

} // ecell4
