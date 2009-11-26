//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-Cell Simulation Environment package
//
//                Copyright (C) 1996-2005 Keio University
//                Copyright (C) 2005 The Molecular Sciences Institute
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// E-Cell is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// E-Cell is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with E-Cell -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Koichi Takahashi <shafi@e-cell.org>,
// E-Cell Project.
//

#ifndef __EVENTSCHEDULER_HPP
#define __EVENTSCHEDULER_HPP

#include "DynamicPriorityQueue.hpp"

namespace libecs
{

    /** @file */
    
    /**
       EventBase
       
       void fire()
       {
         // do what this event is supposed to do.
       }

    */

    class EventBase
    {
        
    public:
        
        EventBase( const double time )
            :
            time( time )
        {
            ; // do nothing
        }

        void setTime( const double time )
        {
            this->time = time;
        }

        const double getTime() const
        {
            return this->time;
        }
   

        const bool operator<= ( const EventBase& rhs ) const
        {
            if( getTime() <= rhs.getTime() )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        const bool operator< ( const EventBase& rhs ) const
        {
            if( getTime() < rhs.getTime() )
            {
                return true;
            }
            else
            {
                return false;
            }
        }


        const bool operator== ( const EventBase& rhs ) const
        {
            if( getTime() == rhs.getTime() )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        const bool operator!= ( const EventBase& rhs ) const
        {
            return ! this->operator==( rhs );
        }



        // dummy, because DynamicPriorityQueue requires this. better without.
        EventBase()
            :
            time( -1.0 )
        {
            ; // do nothing
        }


    private:

        double             time;
    };


    /**
       Event scheduler.

       This class works as a sequential
       event scheduler with a heap-tree based priority queue.

    */

    template <class Event_>
    class EventScheduler
    {
          
    public:

        typedef Event_ Event;
        typedef DynamicPriorityQueue<Event> EventPriorityQueue;

        typedef typename DynamicPriorityQueue<Event>::Index EventIndex;
        typedef typename DynamicPriorityQueue<Event>::ID EventID;

//      typedef std::vector<EventIndex> EventIndexVector;
//      typedef std::vector<EventIndexVector> EventIndexVectorVector;


        EventScheduler()
            :
            time( 0.0 )
        {
            ; // do nothing
        }

        ~EventScheduler()
        {
            ; // do nothing
        }


        const double getTime() const
        {
            return time;
        }

        const double getTopTime() const
        {
            assert( this->getSize() != 0 );  // FIXME: use exception

            return getTopEvent().getTime();
        }

        const EventIndex getSize() const
        {
            return this->eventPriorityQueue.getSize();
        }

        const Event& getTopEvent() const
        {
            return this->eventPriorityQueue.getTop();
        }

        const Event popTopEvent()
        {
            const Event topEvent(getTopEvent());
            this->eventPriorityQueue.popTop();
            return topEvent;
        }

        EventID getTopID() const
        {
            return this->eventPriorityQueue.getTopID();
        }

        const Event& peekSecondEvent() const
        {
            return this->eventPriorityQueue.peekSecond();
        }

        const Event& getEvent( const EventID id ) const
        {
            return this->eventPriorityQueue.get( id );
        }

        const Event& getEventByIndex( const EventIndex index ) const
        {
            return this->eventPriorityQueue.getByIndex( index );
        }

        void step()
        {
            Event topEvent( popTopEvent() );
            this->time = topEvent.getTime();

            topEvent.fire();
        }

        void clear()
        {
            time = 0.0;
            this->eventPriorityQueue.clear();
        }

        const EventID addEvent( const Event& event )
        {
            return this->eventPriorityQueue.push( event );
        }

        void removeEvent( const EventID id )
        {
            this->eventPriorityQueue.pop( id );
        }


        void updateEventTime( const EventID id, const double t )
        {
            const EventIndex index( this->eventPriorityQueue.getIndex( id ) );
            Event& event( this->eventPriorityQueue.getByIndex( index ) );

            event.setTime( t );
            this->eventPriorityQueue.move( index );
        }

        const bool check() const
        {
            return this->eventPriorityQueue.check();
        }

    private:

        EventPriorityQueue       eventPriorityQueue;

        double                   time;


    };

    /*@}*/

} // namespace libecs




#endif /* __EVENTSCHEDULER_HPP */




/*
  Do not modify
  $Author: shafi $
  $Revision: 2529 $
  $Date: 2005-11-19 01:36:40 -0800 (Sat, 19 Nov 2005) $
  $Locker$
*/

