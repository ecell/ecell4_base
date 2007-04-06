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
    
    def testInstantiation(self):
        scheduler = mod.EventScheduler()
        self.failIf( scheduler == None )

    def testEmptyState(self):
        scheduler = mod.EventScheduler()
        self.failIf( scheduler.getSize() != 0 )
        self.failIf( scheduler.getTime() != 0.0 )
        # what if getTopEvent() are called here?

    def testOneEvent(self):
        scheduler = mod.EventScheduler()

        event = TestEvent()
        scheduler.addEvent( 0.0, event )
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getNextTime() != 0.0 )
        self.failIf( scheduler.getTopEvent()[1] != event )
        
        scheduler.step()
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getNextTime() != 1.0 )
        self.failIf( scheduler.getNextTime() != scheduler.getTopEvent()[0] )
        self.failIf( scheduler.getTopEvent()[1] != event )

    def testEventPop(self):
        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event1.dt = - 1.0

        scheduler.addEvent( 0.0, event1 )
        self.failIf( scheduler.getSize() != 1 )
        
        scheduler.step()

        self.failIf( scheduler.getSize() != 0 )
        self.failIf( scheduler.getTime() != 0.0 )

        
    def testEventPop2(self):

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
        self.failIf( scheduler.getTopEvent()[1] != event1 )
        self.failIf( scheduler.getNextTime() != 1.0 )


    # not supported
    def testEventCreationDuringStepping(self):

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
        self.failIf( scheduler.getTopEvent()[1] == event1 or
                     scheduler.getTopEvent()[1] == event2 )
        self.failIf( scheduler.getNextTime() != 0.5 )


    # not supported
    def testEventCreationDuringSteppingWithPop(self):

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
        self.failIf( scheduler.getTopEvent()[1] == event1 or
                     scheduler.getTopEvent()[1] == event2 )
        self.failIf( scheduler.getNextTime() != 0.5 )


if __name__ == "__main__":
    unittest.main()
