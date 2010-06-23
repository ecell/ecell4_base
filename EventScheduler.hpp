//
// written by Koichi Takahashi <shafi@e-cell.org>,
// E-Cell Project.
//

#ifndef __EVENTSCHEDULER_HPP
#define __EVENTSCHEDULER_HPP

#include "DynamicPriorityQueue.hpp"

/**
   Event scheduler.

   This class works as a sequential
   event scheduler with a heap-tree based priority queue.

*/

template <class Tevent_>
class EventScheduler
{
public:
    typedef Tevent_ event_type;
    typedef typename event_type::time_type time_type;

protected:
    struct event_comparator
    {
        bool operator()(event_type const& lhs, event_type const& rhs) const
        {
            return lhs.time() <= rhs.time();
        }
    };

    typedef DynamicPriorityQueue<event_type, event_comparator> EventPriorityQueue;

public:
    typedef typename EventPriorityQueue::size_type size_type;
    typedef typename EventPriorityQueue::identifier_type identifier_type;
    typedef typename EventPriorityQueue::value_type value_type;

public:


    EventScheduler()
        : time_( 0.0 ) {}

    ~EventScheduler() {}

    time_type time() const
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
        const value_type top(eventPriorityQueue_.top());
        eventPriorityQueue_.pop();
        return top;
    }

    value_type const& second() const
    {
        return eventPriorityQueue_.second();
    }

    event_type const& get(identifier_type const& id) const
    {
        return eventPriorityQueue_.get(id);
    }

    void step()
    {
        value_type topEvent(pop());
        time_ = topEvent.second.time();
        topEvent.second.fire();
    }

    void clear()
    {
        time_ = 0.0;
        eventPriorityQueue_.clear();
    }

    identifier_type add(event_type const& event)
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

private:
    EventPriorityQueue eventPriorityQueue_;
    time_type time_;
};

#endif /* __EVENTSCHEDULER_HPP */
