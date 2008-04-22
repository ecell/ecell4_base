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
       
       A subclass must define three customization points;
    
       void fire()
       {
       (1) do what this event is supposed to do.
       (2) setTime( next scheduled time of this event );
       }

       void update( const Event& anEvent )
       {
       Given the last fired Event (anEvent) that this Event
       depends on,

       (1) recalculate scheduled time (if necessary).
       (2) setTime( new scheduled time ).
       }

       const bool isDependentOn( const Event& anEvent )
       {
       Return true if this Event must be updated when the
       given Event (anEvent) fired.  Otherwise return false;
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

//	typedef std::vector<EventIndex> EventIndexVector;
//	typedef std::vector<EventIndexVector> EventIndexVectorVector;


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

	EventID getTopID() const
	{
	    return this->eventPriorityQueue.getTopID();
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

	    // Here I copy construct the top event and use its event
	    // ID to reschedule it.  This is necessary if events can
	    // be created or deleted within fire() and the dynamic
	    // priority queue can reallocate internal data structures.
	    // Most of the cost of using this is optimized away when
	    // the dynamic priority queue has a VolatileIDPolicy.
	    Event topEvent( getTopEvent() );
	    const EventID ID( this->eventPriorityQueue.getTopID() );
	    this->time = topEvent.getTime();

	    // Fire top
	    topEvent.fire();

	    // If the event is rescheduled into the past, remove it.
	    // Otherwise, reuse the event.
	    if( topEvent.getTime() >= getTime() )
	    {
		this->eventPriorityQueue.replace( ID, topEvent );
	    }
	    else
	    {
		this->eventPriorityQueue.pop( ID );
	    }

//	    assert( getNextTime() >= getTime() );

	    // update dependent events
//	    const EventIndexVector&
//		anEventIndexVector( this->eventDependencyArray[ topEventIndex ] );

/*
	    for( typename EventIndexVector::const_iterator 
		     i( anEventIndexVector.begin() );
		 i != anEventIndexVector.end(); ++i )
	    {
		const EventIndex anIndex( *i );

		updateEvent( anIndex, currentTime );
	    }
*/
	}

/*
	void updateAllEvents( const double aCurrentTime )
	{
	    const EventIndex aSize( getSize() );
	    for( EventIndex anIndex( 0 ); anIndex != aSize; ++anIndex )
	    {
		updateEvent( anIndex, aCurrentTime );
	    }
	}

	void updateEvent( const EventIndex anIndex, const double aCurrentTime )
	{
	    Event& anEvent( this->eventPriorityQueue.getIndex( anIndex ) );
	    const double anOldTime( anEvent.getTime() );
	    anEvent.update( aCurrentTime );
	    const double aNewTime( anEvent.getTime() );

	    // this->eventPriorityQueue.move( anIndex );
	    if( aNewTime >= anOldTime )
	    {
		this->eventPriorityQueue.moveDown( anIndex );
	    }
	    else
	    {
		this->eventPriorityQueue.moveUp( anIndex );
	    }
	}
*/

//	void updateAllEventDependency();  // update all

//	void updateEventDependency( const EventIndex anIndex );
    
	void clear()
	{
	    this->eventPriorityQueue.clear();
//	    this->eventDependencyArray.clear();
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



	// this is here for DiscreteEventStepper::log().
	// should be removed in future. 
//	const EventIndexVector& getDependencyVector( const EventIndex anIndex )
//	{
//	    return this->eventDependencyArray[ anIndex ] ;
//	}

        const bool check() const
        {
            return this->eventPriorityQueue.check();
        }

    private:

	EventPriorityQueue       eventPriorityQueue;

//	EventIndexVectorVector   eventDependencyArray;

	double                   time;


    };

  
/*

    template < class Event >
    void EventScheduler<Event>::updateAllEventDependency()
    {
	this->eventDependencyArray.resize( this->eventPriorityQueue.getSize() );
    
	for( EventIndex i1( 0 ); i1 != this->eventPriorityQueue.getSize(); ++i1 )
	{
	    updateEventDependency( i1 );
	}
    }

    template < class Event >
    void EventScheduler<Event>::
    updateEventDependency( const EventIndex i1 )
    {
	const Event& anEvent1( this->eventPriorityQueue[ i1 ] );

	EventIndexVector& anEventIndexVector( this->eventDependencyArray[ i1 ] );
	anEventIndexVector.clear();

	for( EventIndex i2( 0 ); i2 < this->eventPriorityQueue.getSize(); ++i2 )
	{
	    if( i1 == i2 )
	    {
		// don't include itself
		continue;
	    }
	
	    const Event& anEvent2( this->eventPriorityQueue[ i2 ] );
	
	    if( anEvent2.isDependentOn( anEvent1 ) )
	    {
		anEventIndexVector.push_back( i2 );
	    }
	}
    
	std::sort( anEventIndexVector.begin(), anEventIndexVector.end() );
    }
*/




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

