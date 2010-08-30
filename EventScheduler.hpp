#ifndef __EVENTSCHEDULER_HPP
#define __EVENTSCHEDULER_HPP
//
// written by Koichi Takahashi <shafi@e-cell.org>,
// E-Cell Project.
//


#include "DynamicPriorityQueue.hpp"
#include <boost/shared_ptr.hpp>
#include <stdexcept>

/**
   Event scheduler.

   This class works as a sequential
   event scheduler with a heap-tree based priority queue.

*/

template<typename Ttime_>
class EventScheduler
{
public:
    typedef Ttime_ time_type;

    struct Event
    {
        typedef Ttime_ time_type;

        Event(time_type const& time): time_(time) {}

        virtual ~Event() {}

        time_type const& time() const
        {
            return time_;
        }

    protected:
        const time_type time_;
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

    typedef DynamicPriorityQueue<boost::shared_ptr<Event>, event_comparator> EventPriorityQueue;

public:
    typedef typename EventPriorityQueue::size_type size_type;
    typedef typename EventPriorityQueue::identifier_type identifier_type;
    typedef typename EventPriorityQueue::value_type value_type;
    typedef boost::iterator_range<typename EventPriorityQueue::const_iterator> events_range;

public:


    EventScheduler(): time_( 0.0 ) {}

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
        return boost::make_iterator_range(eventPriorityQueue_.begin(),
                                          eventPriorityQueue_.end());
    }

private:
    EventPriorityQueue eventPriorityQueue_;
    time_type time_;
};

#endif /* __EVENTSCHEDULER_HPP */
