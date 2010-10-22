#!/usr/bin/env python

import unittest

import weakref

import _gfrd as mod


class Delegate(object):

    def __init__(self, obj, method):
        self.obj = weakref.proxy(obj)
        self.method = method

    def __call__(self, arg):
        return self.method(self.obj, arg)



class EventSchedulerTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass
    
    def test_instantiation(self):
        scheduler = mod.EventScheduler()
        self.failIf(scheduler == None)

    def test_empty_state(self):
        scheduler = mod.EventScheduler()
        self.failIf(scheduler.size != 0)
        self.failIf(scheduler.time != 0.0)
        # what if getTopEvent() are called here?

    def test_one_event(self):
        scheduler = mod.EventScheduler()

        event = mod.Event(1.0)
        id = scheduler.add(event)
        self.failIf(scheduler.time != 0.0)
        self.failIf(scheduler.top[1].time != 1.0)
        self.failIf(scheduler.top[0] != id)

        self.assertEqual((id, event), scheduler.pop())
        self.failIf(scheduler.size != 0)
        self.failIf(scheduler.time != 1.0)


    def test_peek_second_event(self):

        scheduler = mod.EventScheduler()

        event1 = mod.Event(1.0)
        event2 = mod.Event(0.5)

        event1_id = scheduler.add(event1)

        id2 = scheduler.add(event2)

        self.assertEqual(2, scheduler.size)

        second = scheduler.second

        self.assertEqual(1.0, second[1].time)
        self.assertEqual(event1_id, second[0])




if __name__ == "__main__":
    unittest.main()
