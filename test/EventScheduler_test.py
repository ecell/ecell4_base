#!/usr/bin/env python

import unittest

import _gfrd as mod



class TestEvent:

    def __init__( self ):
        self.fired = 0
        self.dt = 1.0

    def fire( self ):
        self.fired += 1
        return self.dt


class TestEvent2:

    def __init__( self, scheduler ):
        self.fired = 0
        self.dt = 1.0
        self.scheduler = scheduler

    def fire( self ):
        self.fired += 1
        self.scheduler.addEvent( self.scheduler.getTime(), TestEvent() )
        return self.dt

    

class EventSchedulerTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test1Instantiation(self):
        scheduler = mod.EventScheduler()
        self.failIf( scheduler == None )

    def test2EmptyState(self):
        scheduler = mod.EventScheduler()
        self.failIf( scheduler.getSize() != 0 )
        self.failIf( scheduler.getTime() != 0.0 )
        # what if getTopEvent() are called here?

    def test3OneEvent(self):
        scheduler = mod.EventScheduler()

        event = TestEvent()
        scheduler.addEvent( 0.0, event )
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getNextTime() != 0.0 )
        self.failIf( scheduler.getTopEvent().getObj() != event )
        
        scheduler.step()
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getNextTime() != 1.0 )
        self.failIf( scheduler.getNextTime() != 
                     scheduler.getTopEvent().getTime() )
        self.failIf( scheduler.getTopEvent().getObj() != event )

    def test4EventPop(self):
        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event1.dt = - 1.0

        scheduler.addEvent( 0.0, event1 )
        self.failIf( scheduler.getSize() != 1 )
        
        scheduler.step()

        self.failIf( scheduler.getSize() != 0 )
        self.failIf( scheduler.getTime() != 0.0 )

        
    def test5EventPop2(self):

        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event2 = TestEvent()
        event1.dt = 1.0
        event2.dt = -1.0

        scheduler.addEvent( 0.0, event1 )
        scheduler.addEvent( 0.5, event2 )
        self.failIf( scheduler.getSize() != 2 )
        
        scheduler.step()
        self.failIf( scheduler.getTime() != 0.0 )
        scheduler.step()

        self.failIf( scheduler.getTime() != 0.5 )
        self.failIf( scheduler.getSize() != 1 )
        self.failIf( scheduler.getTopEvent().getObj() != event1 )
        self.failIf( scheduler.getNextTime() != 1.0 )


    def test6EventCreationDuringStepping(self):

        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event2 = TestEvent2( scheduler )
        event1.dt = 1.0
        event2.dt = 1.0

        scheduler.addEvent( 0.0, event1 )
        scheduler.addEvent( 0.5, event2 )
        self.failIf( scheduler.getSize() != 2 )
        
        scheduler.step()
        self.failIf( scheduler.getTime() != 0.0 )
        scheduler.step()

        self.failIf( scheduler.getTime() != 0.5 )
        self.failIf( scheduler.getSize() != 3 )
        self.failIf( scheduler.getTopEvent().getObj() == event1 or
                     scheduler.getTopEvent().getObj() == event2 )
        self.failIf( scheduler.getNextTime() != 0.5 )


    def test7EventCreationDuringSteppingWithPop(self):

        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event2 = TestEvent2( scheduler )
        event1.dt = 1.0
        event2.dt = - 1.0

        scheduler.addEvent( 0.0, event1 )
        scheduler.addEvent( 0.5, event2 )
        self.failIf( scheduler.getSize() != 2 )
        
        scheduler.step()
        self.failIf( scheduler.getTime() != 0.0 )
        scheduler.step()

        self.failIf( scheduler.getTime() != 0.5 )
        self.failIf( scheduler.getSize() != 2 )
        self.failIf( scheduler.getTopEvent().getObj() == event1 or
                     scheduler.getTopEvent().getObj() == event2 )
        self.failIf( scheduler.getNextTime() != 0.5 )


if __name__ == "__main__":
    unittest.main()
