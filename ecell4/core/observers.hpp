#ifndef __ECELL4_OBSERVER_HPP
#define __ECELL4_OBSERVER_HPP

#include "types.hpp"
#include "functions.hpp"
#include "Space.hpp"
#include "Simulator.hpp"

#include <fstream>
#include <boost/format.hpp>
#include <time.h>


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

    virtual const Real next_time() const;
    virtual void initialize(const boost::shared_ptr<Space>& space);
    virtual void finalize(const boost::shared_ptr<Space>& space);
    virtual void reset();

    virtual bool fire(const Simulator* sim, const boost::shared_ptr<Space>& space);
    // virtual bool fire(const Simulator* sim, const boost::shared_ptr<Space>& space) = 0;

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

    const Real next_time() const;
    const Integer num_steps() const;
    const Integer count() const;
    virtual void initialize(const boost::shared_ptr<Space>& space);
    virtual bool fire(const Simulator* sim, const boost::shared_ptr<Space>& space);
    virtual void reset();

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

    void log(const boost::shared_ptr<Space>& space);

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

    virtual void initialize(const boost::shared_ptr<Space>& space);
    virtual bool fire(const Simulator* sim, const boost::shared_ptr<Space>& space);
    virtual void reset();
    NumberLogger::data_container_type data() const;
    NumberLogger::species_container_type targets() const;

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
        : base_type(true), logger_(species), num_steps_(0)
    {
        ;
    }

    virtual ~NumberObserver()
    {
        ;
    }

    virtual void initialize(const boost::shared_ptr<Space>& space);
    virtual void finalize(const boost::shared_ptr<Space>& space);
    virtual bool fire(const Simulator* sim, const boost::shared_ptr<Space>& space);
    virtual void reset();
    const Integer num_steps() const;
    NumberLogger::data_container_type data() const;
    NumberLogger::species_container_type targets() const;

protected:

    NumberLogger logger_;
    Integer num_steps_;
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

    const Real next_time() const;

    const Integer num_steps() const
    {
        return num_steps_;
    }

    virtual void initialize(const boost::shared_ptr<Space>& space);
    virtual bool fire(const Simulator* sim, const boost::shared_ptr<Space>& space);
    virtual void reset();

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

    virtual void initialize(const boost::shared_ptr<Space>& space);
    virtual bool fire(const Simulator* sim, const boost::shared_ptr<Space>& space);
    virtual void reset();
    NumberLogger::data_container_type data() const;
    NumberLogger::species_container_type targets() const;

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

    virtual void initialize(const boost::shared_ptr<Space>& space);
    virtual bool fire(const Simulator* sim, const boost::shared_ptr<Space>& space);

    inline const std::string filename() const
    {
        return filename(num_steps());
    }

    const std::string filename(const Integer idx) const;

    const std::string& prefix() const
    {
        return prefix_;
    }

protected:

    std::string prefix_;
};

class FixedIntervalCSVObserver
    : public FixedIntervalObserver
{
public:

    typedef FixedIntervalObserver base_type;

    typedef std::vector<std::pair<ParticleID, Particle> >
        particle_container_type;
    typedef utils::get_mapper_mf<Species::serial_type, unsigned int>::type
        serial_map_type;

public:

    FixedIntervalCSVObserver(
        const Real& dt, const std::string& filename)
        : base_type(dt), prefix_(filename), species_()
    {
        ;
    }

    FixedIntervalCSVObserver(
        const Real& dt, const std::string& filename,
        const std::vector<std::string>& species)
        : base_type(dt), prefix_(filename), species_(species)
    {
        ;
    }

    virtual ~FixedIntervalCSVObserver()
    {
        ;
    }

    virtual void initialize(const boost::shared_ptr<Space>& space);
    virtual bool fire(const Simulator* sim, const boost::shared_ptr<Space>& space);
    void write_particles(
        std::ofstream& ofs, const particle_container_type& particles,
        const Species::serial_type label = "");
    void log(const boost::shared_ptr<Space>& space);
    const std::string filename() const;
    virtual void reset();

protected:

    std::string prefix_;
    std::vector<std::string> species_;
    serial_map_type serials_;
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
        : base_type(dt), pids_(pids), resolve_boundary_(resolve_boundary),
        trajectories_(pids.size()), strides_(pids.size()), t_()
    {
        ;
    }

    FixedIntervalTrajectoryObserver(
        const Real& dt, const bool resolve_boundary = true)
        : base_type(dt), pids_(), resolve_boundary_(resolve_boundary),
        trajectories_(), strides_(), t_()
    {
        ;
    }

    virtual ~FixedIntervalTrajectoryObserver()
    {
        ;
    }

    virtual void initialize(const boost::shared_ptr<Space>& space);
    virtual bool fire(const Simulator* sim, const boost::shared_ptr<Space>& space);
    virtual void reset();

    const std::vector<std::vector<Real3> >& data() const;
    const Integer num_tracers() const;
    const std::vector<Real>& t() const;

protected:

    std::vector<ParticleID> pids_;
    bool resolve_boundary_;
    std::vector<std::vector<Real3> > trajectories_;
    std::vector<Real3> strides_;
    std::vector<Real> t_;
};


class TimeoutObserver
    : public Observer
{
public:

    typedef Observer base_type;

public:

    TimeoutObserver(const Real interval)
        : base_type(true), interval_(interval), duration_(0.0)
    {
        ;
    }

    TimeoutObserver()
        : base_type(true), interval_(inf), duration_(0.0), acc_(0.0)
    {
        ;
    }

    virtual ~TimeoutObserver()
    {
        ;
    }

    virtual void initialize(const boost::shared_ptr<Space>& space);
    virtual void finalize(const boost::shared_ptr<Space>& space);
    virtual bool fire(const Simulator* sim, const boost::shared_ptr<Space>& space);
    virtual void reset();

    const Real interval() const
    {
        return interval_;
    }

    const Real duration() const
    {
        return duration_;
    }

    const Real accumulation() const
    {
        return acc_;
    }

protected:

    Real interval_;
    Real duration_;
    Real acc_;
    time_t tstart_;
};

} // ecell4

#endif /* __ECELL4_OBSEVER_HPP */
