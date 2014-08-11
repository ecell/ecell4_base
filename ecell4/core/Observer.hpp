#ifndef __ECELL4_OBSERVER_HPP
#define __ECELL4_OBSERVER_HPP

#include "types.hpp"
#include "Space.hpp"


namespace ecell4
{

class Observer
{
public:

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

    virtual void fire(const Space* space) = 0;
};

class FixedIntervalObserver
    : public Observer
{
public:

    FixedIntervalObserver(const Real& dt)
        : tnext_(0.0), dt_(dt), num_steps_(0)
    {
        ;
    }

    virtual ~FixedIntervalObserver()
    {
        ;
    }

    const Real next_time() const
    {
        return tnext_;
    }

    virtual void initialize(const Space* space)
    {
        tnext_ = space->t();
        num_steps_ = 0;
    }

    virtual void fire(const Space* space)
    {
        tnext_ += dt_;
        ++num_steps_;
    }

protected:

    Real tnext_, dt_;
    Integer num_steps_;
};

class FixedIntervalNumberObserver
    : public FixedIntervalObserver
{
public:

    typedef FixedIntervalObserver base_type;
    typedef std::vector<std::vector<Real> > data_container_type;
    typedef std::vector<Species> species_container_type;

public:

    FixedIntervalNumberObserver(const Real& dt, const std::vector<std::string> species)
        : base_type(dt)
    {
        species_.reserve(species.size());
        for (std::vector<std::string>::const_iterator i(species.begin());
            i != species.end(); ++i)
        {
            species_.push_back(Species(*i));
        }
    }

    virtual ~FixedIntervalNumberObserver()
    {
        ;
    }

    virtual void initialize(const Space* space)
    {
        base_type::initialize(space);
        data_.clear();
    }

    virtual void fire(const Space* space)
    {
        data_container_type::value_type tmp;
        tmp.push_back(space->t());
        for (species_container_type::const_iterator i(species_.begin());
            i != species_.end(); ++i)
        {
            tmp.push_back(space->num_molecules(*i));
        }
        data_.push_back(tmp);

        base_type::fire(space);
    }

    data_container_type data() const
    {
        return data_;
    }

protected:

    data_container_type data_;
    species_container_type species_;
};

} // ecell4

#endif /* __ECELL4_OBSEVER_HPP */
