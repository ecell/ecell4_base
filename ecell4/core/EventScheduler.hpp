#ifndef __ECELL4_EVENTSCHEDULER_HPP
#define __ECELL4_EVENTSCHEDULER_HPP

#include <boost/range/iterator_range.hpp>
#include <boost/shared_ptr.hpp>
#include <stdexcept>

#include "types.hpp"
#include "DynamicPriorityQueue.hpp"


namespace ecell4
{

class EventScheduler
{
public:

    struct Event
    {
    public:

        Event(Real const& time) : time_(time) {}

        virtual ~Event() {}

        virtual void fire() {}

        Real const& time() const
        {
            return time_;
        }

        //XXX: deprecate me
        Real const& dt() const
        {
            return dt_;
        }

        virtual void interrupt(Real const& t) {}

    protected:

        Real time_;
        //XXX: deprecate me
        Real dt_;
    };

protected:

    struct event_comparator
    {
        bool operator()(boost::shared_ptr<Event> const& lhs,
                boost::shared_ptr<Event> const& rhs) const
        {
            return lhs->time() <= rhs->time();
        }
    };

    typedef DynamicPriorityQueue<boost::shared_ptr<Event>, event_comparator>
        EventPriorityQueue;

public:

    typedef EventPriorityQueue::size_type size_type;
    typedef EventPriorityQueue::identifier_type identifier_type;
    typedef EventPriorityQueue::value_type value_type;
    typedef boost::iterator_range<EventPriorityQueue::const_iterator>
        events_range;

public:

    EventScheduler() : time_(0.0) {}

    ~EventScheduler() {}

    Real time() const
    {
        return time_;
    }

    size_type size() const
    {
        return eventPriorityQueue_.size();
    }

    value_type const& top() const
    {
        return eventPriorityQueue_.top();
    }

    value_type pop()
    {
        if (eventPriorityQueue_.empty())
        {
            throw std::out_of_range("queue is empty");
        }
        const value_type top(eventPriorityQueue_.top());
        eventPriorityQueue_.pop();
        time_ = top.second->time();
        return top;
    }

    value_type const& second() const
    {
        return eventPriorityQueue_.second();
    }

    boost::shared_ptr<Event> get(identifier_type const& id) const
    {
        return eventPriorityQueue_.get(id);
    }

    void clear()
    {
        time_ = 0.0;
        eventPriorityQueue_.clear();
    }

    identifier_type add(boost::shared_ptr<Event> const& event)
    {
        return eventPriorityQueue_.push(event);
    }

    void remove(identifier_type const& id)
    {
        eventPriorityQueue_.pop(id);
    }

    void update(value_type const& pair)
    {
        eventPriorityQueue_.replace(pair);
    }

    bool check() const
    {
        return eventPriorityQueue_.check();
    }

    events_range events() const
    {
        return boost::make_iterator_range(
            eventPriorityQueue_.begin(), eventPriorityQueue_.end());
    }

    const Real next_time() const
    {
        if (size() > 0)
        {
            return top().second->time();
        }
        else
        {
            return inf;
        }
    }

protected:

    EventPriorityQueue eventPriorityQueue_;
    Real time_;
};

} // ecell4

#endif /* __ECELL4_EVENTSCHEDULER_HPP */
