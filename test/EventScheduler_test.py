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

class TestEvent3:

    def __init__( self, scheduler, peer ):
        self.fired = 0
        self.dt = 1.0
        self.scheduler = scheduler
        self.peer = peer

    def fire( self ):
        self.fired += 1
        self.scheduler.updateEvent( self.peer.id, self.scheduler.getTime(),
                                    self.peer )
        return self.dt

class TestEvent4:

    def __init__( self, scheduler ):
        self.fired = 0
        self.dt = 1.0
        self.scheduler = scheduler

    def fire( self ):
        self.fired += 1
        self.scheduler.addEvent( self.scheduler.getTime(), TestEvent() )
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
        self.failIf( scheduler.getTopEvent().getObj() != event )
        
        scheduler.step()
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getNextTime() != 1.0 )
        self.failIf( scheduler.getNextTime() != 
                     scheduler.getTopEvent().getTime() )
        self.failIf( scheduler.getTopEvent().getObj() != event )

    def testTwoEventsSameTime1(self):
        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event2 = TestEvent()
        event1.dt = 1.0
        event2.dt = 0.5
        scheduler.addEvent( 0.0, event1 )
        scheduler.addEvent( 0.0, event2 )
        self.assertEqual( 2, scheduler.getSize() )
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getNextTime() != 0.0 )
        
        scheduler.step()
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getNextTime() != 0.0 )

        scheduler.step()
        self.failIf( scheduler.getNextTime() != 0.5 )
        self.failIf( scheduler.getNextTime() != 
                     scheduler.getTopEvent().getTime() )
        self.failIf( scheduler.getTopEvent().getObj() != event2 )


    def testTwoEventsSameTime2(self):
        scheduler = mod.EventScheduler()

        event1 = TestEvent4( scheduler )
        event1.dt = 0.5
        scheduler.addEvent( 0.0, event1 )
        self.assertEqual( 1, scheduler.getSize() )
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getNextTime() != 0.0 )
        
        scheduler.step()
        self.assertEqual( 3, scheduler.getSize() )
        self.failIf( scheduler.getTime() != 0.0 )
        self.failIf( scheduler.getNextTime() != 0.0 )

        scheduler.step()
        scheduler.step()
        self.assertEqual( 0.5, scheduler.getNextTime() )
        self.assertEqual( event1, scheduler.getTopEvent().getObj() )



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
        self.failIf( scheduler.getTopEvent().getObj() != event1 )
        self.failIf( scheduler.getNextTime() != 1.0 )


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
        self.failIf( scheduler.getTopEvent().getObj() == event1 or
                     scheduler.getTopEvent().getObj() == event2 )
        self.failIf( scheduler.getNextTime() != 0.5 )


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
        self.failIf( scheduler.getTopEvent().getObj() == event1 or
                     scheduler.getTopEvent().getObj() == event2 )
        self.failIf( scheduler.getNextTime() != 0.5 )


    def testEventUpdateDuringStepping(self):

        scheduler = mod.EventScheduler()

        event1 = TestEvent()
        event2 = TestEvent3( scheduler, event1 ) 
        event1.dt = 1.0
        event2.dt = 1.0

        event1.id = scheduler.addEvent( 0.0, event1 )
        scheduler.addEvent( 0.5, event2 )
        self.assertEqual( 2, scheduler.getSize() )
        
        self.failIf( scheduler.getTopEvent().getObj() != event1 )
        scheduler.step()
        self.assertEqual( 0.0, scheduler.getTime() )

        # event2 updates event1 to step immediately
        self.failIf( scheduler.getTopEvent().getObj() != event2 )
        scheduler.step()

        self.assertEqual( 0.5, scheduler.getTime() )
        self.assertEqual( 2, scheduler.getSize() )
        self.failIf( scheduler.getTopEvent().getObj() != event1 )
        self.failIf( scheduler.getNextTime() != 0.5 )




if __name__ == "__main__":
    unittest.main()
