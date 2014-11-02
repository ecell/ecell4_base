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
        : base_type(false), t0_(0.0), dt_(dt), num_steps_(0)
        // : base_type(false), tnext_(0.0), dt_(dt), num_steps_(0)
    {
        ;
    }

    virtual ~FixedIntervalObserver()
    {
        ;
    }

    const Real next_time() const
    {
        return t0_ + dt_ * num_steps_;
        // return tnext_;
    }

    const Integer num_steps() const
    {
        return num_steps_;
    }

    virtual void initialize(const Space* space)
    {
        t0_ = space->t();
        // tnext_ = space->t();
        num_steps_ = 0;
    }

    virtual void fire(const Simulator* sim, const Space* space)
    {
        // tnext_ += dt_;
        ++num_steps_;
    }

protected:

    Real t0_, dt_;
    // Real tnext_, dt_;
    Integer num_steps_;
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

} // ecell4

#endif /* __ECELL4_OBSEVER_HPP */
