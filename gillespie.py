#!/usr/env python

import numpy
import weakref
import logging
import os

from _gfrd import * # FIX ME
from constants import EventType # FIX ME
import utils
import logger
import myrandom

log = logging.getLogger('gillespie')

def rng_uniform():
    '''Returns a positive random number
    '''
    while True:
        rng = myrandom.uniform(0.0, 1.0)
        if rng > 0.0: return rng

class Logger(logger.Logger):

    def __init__(self, logname='log', directory='data', comment=''):
        logger.Logger.__init__(self, logname, directory, comment)

    def prepare_timecourse_file(self, simulator):
        if not os.path.exists(self.directory):
            os.mkdir(self.directory)
        timecourse_filename = '%s_tc.dat' % self.logname
        self.timecourse_file = open(
            os.path.join(self.directory, timecourse_filename), 'w')
        self.write_timecourse_comment(self.comment)

        species_name_list = '\'' + \
            "\', \'".join(str(id) for id in simulator.get_species_id()) + '\''
        columns = '[\'t\', ' + species_name_list + ']'
        self.write_timecourse_comment('@ columns= ' + columns)

    def write_timecourse(self, simulator):
        data = []
        self.timecourse_file.write('%g\t' % simulator.t)
        self.timecourse_file.write('\t'.join(
                str(simulator.get_pool_size(id))
                for id in simulator.get_species_id()) + '\n')
        self.timecourse_file.flush()

    def write_particles(self, simulator):
        # dummy
        pass

class Delegate(object):

    def __init__(self, obj, method):
        self.ref = weakref.ref(obj)
        self.method = method

    def __call__(self, *arg):
        return self.method(self.ref(), *arg)

class GillespieEvent(Event):
    __slot__ = ['func', 'rr']
    def __init__(self, time, func, rr):
        Event.__init__(self, time)
        self.func = func
        self.rr = rr

class ReactionRuleCache(object):

    def __init__(self, rr, reactants, products, k):
        self.rr = rr

        self.reactants = reactants
        self.products = products
        self.k = k

        self.eventID = None
        self.eventType = None

    def accessors(self):
        if self.eventType == EventType.SINGLE_REACTION: # FIX ME
            return self.reactants
        else:
            return []

    def mutators(self):
        if self.eventType == EventType.SINGLE_REACTION: # FIX ME
            return self.reactants + self.products
        else:
            return []

    def is_dependent_on(self, rr):
        for id1 in self.accessors():
            if id1 in rr.mutators():
                return True

        return False

class GillespieSimulatorBase(object):

    def __init__(self):
        self.model = None
        self.speciesDict = {}
        self.stateArray = numpy.array([])

        self.set_volume(utils.INF)

    def set_model(self, model):
        self.model = model
        self.network_rules = NetworkRulesWrapper(model.network_rules)

    def initialize(self):
        pass

    def reset(self):
        pass

    def set_world_size(self, size):
        if size == utils.INF:
            self.set_volume(utils.INF)
        else:
            volume = size * size * size
            self.set_volume(volume)

    def get_world_size(self):
        return self.volume ** (1.0 / 3.0)

    def set_volume(self, volume):
        self.volume = volume

    def get_volume(self):
        return self.volume

    def create_reaction_rule_cache(self, rr):
        reactants = [id for id in rr.reactants]
        products = [id for id in rr.products]

        if len(reactants) == 1:
            k = float(rr.k)
        elif len(reactants) == 2:
            st1 = self.model.get_species_type_by_id(reactants[0])
            st2 = self.model.get_species_type_by_id(reactants[1])
            D = float(st1['D']) + float(st2['D'])
            sigma = float(st1['radius']) + float(st2['radius'])
            kD = utils.k_D(D, sigma)
            k = float(rr.k)
            if kD == 0.0:
                k = 0.0
            elif k != 0.0:
                k = utils.k_on(k, kD)

        return ReactionRuleCache(rr, reactants, products, k)

    def get_reaction_rule1(self, st):
        return self.__get_reaction_rule(st)

    def get_reaction_rule2(self, st1, st2):
        return self.__get_reaction_rule(st1, st2)

    def __get_reaction_rule(self, *args):
        gen = self.network_rules.query_reaction_rule(*args)
        if gen == None:
            return []

        retval = []
        for rr in gen:
            if rr.k == 0:
                continue
            retval.append(self.create_reaction_rule_cache(rr))
        return retval

    def clear(self):
        pass

    def get_pool_size(self, id):
        if id in self.speciesDict.keys():
            return self.stateArray[self.speciesDict[id]]
        else:
            return 0.0

    def remove_particles(self, st, n):
        if __debug__:
            log.info('removing in %s %s particles' % (n, st.id))

        if st.id not in self.speciesDict.keys():
            raise RuntimeError, '%s species doesn\'t exist.' % (st.id)

        self.stateArray[self.speciesDict[st.id]] -= n

    def throw_in_particles(self, st, n):
        if __debug__:
            log.info('throwing in %s %s particles' % (n, st.id))

        if not st.id in self.speciesDict.keys():
            raise RuntimeError, '%s species doesn\'t exist.' % (st.id)

        self.stateArray[self.speciesDict[st.id]] += n

    def add_species(self, id):
        i = len(self.stateArray)
        self.stateArray = numpy.resize(self.stateArray, i + 1)
        self.stateArray[i] = 0.0

        self.speciesDict[id] = i

    def check(self):
        pass

    def dump_population(self):
        buf = ''
#         for id, i in self.speciesDict.items():
#             st = self.model.get_species_type_by_id(id)
#             buf += st['id'] + ':' + str(self.stateArray[i]) + '\t'

        for id, i in self.speciesDict.items():
            buf += str(self.stateArray[i]) + ' '

        return buf

    def get_step_interval(self, rr):
        return self.get_propensity_R(rr) * (- numpy.log(rng_uniform()))

    def get_propensity(self, rr):
        if len(rr.reactants) == 1:
            propensity = self.get_propensity_first_order(rr)

        elif len(rr.reactants) == 2:
            if rr.reactants[0] == rr.reactants[1]:
                propensity \
                    = self.get_propensity_second_order_one_substrate(rr)
            else:
                propensity \
                    = self.get_propensity_second_order_two_substrates(rr)
        
        if propensity < 0.0:
            raise RuntimeError, 'Population size <= -1.0'
            return 0.0

        else:
            return propensity

    def get_propensity_R(self, rr):
        propensity = self.get_propensity(rr)
        if propensity > 0.0:
            return 1.0 / propensity
        else:
            return utils.INF

    def get_propensity_first_order(self, rr):
        value = self.get_pool_size(rr.reactants[0])

        if value > 0.0:
            return rr.k * value
        else:
            return 0.0

    def get_propensity_second_order_two_substrates(self, rr):
        value = self.get_pool_size(rr.reactants[0]) \
            * self.get_pool_size(rr.reactants[1])

        if value > 0.0:
            return rr.k * value / self.volume
        else:
            return 0.0

    def get_propensity_second_order_one_substrate(self, rr):
        value = self.get_pool_size(rr.reactants[0])

        if value > 1.0: # there must be two or more molecules
            return rr.k * 0.5 * value * (value - 1.0) / self.volume
        else:
            return 0.0

class GillespieSimulator(GillespieSimulatorBase):

    def __init__(self, model):
        self.scheduler = EventScheduler()
        GillespieSimulatorBase.__init__(self)

        self.dependencies = {}
        self.last_event = None
        self.last_reaction = None

        self.t = 0.0
        self.dt = 0.0

        self.scheduler.clear()

        volume = model.world_size ** 3
        self.set_volume(volume)
        self.set_model(model)

    def initialize(self):
        GillespieSimulatorBase.initialize(self)

        self.scheduler.clear()

        for id1, i in self.speciesDict.items():
            st1 = self.model.get_species_type_by_id(id1)
            rules = self.get_reaction_rule1(st1)
            for rr in rules:
                self.add_update_event(rr)

            for id2, j in self.speciesDict.items():
                if i > j:
                    continue

                st2 = self.model.get_species_type_by_id(id2)
                rules = self.get_reaction_rule2(st1, st2)
                for rr in rules:
                    self.add_update_event(rr)

    def reset(self):
        GillespieSimulatorBase.reset(self)

    def get_next_time(self):
        if self.scheduler.size == 0:
            return self.t

        return self.scheduler.top[1].time

    def govern(self, id):
        return id in self.speciesDict.keys()

    def stop(self, t):
        if __debug__:
            log.info('stop at %g' % t)

        if self.t == t:
            return

        if t >= self.scheduler.getTopEvent().getTime():
            raise RuntimeError, 'Stop time >= next event time.'

        if t < self.t:
            raise RuntimeError, 'Stop time >= next event time.'

        self.t = t

    def step(self):
        if self.scheduler.size == 0:
            self.t = utils.INF
            self.dt = utils.INF
            self.last_event = None
            self.last_reaction = None
            return

        id, event = self.scheduler.pop()
        self.t, self.last_event = event.time, event

        if self.last_event.rr.eventType == EventType.SINGLE_REACTION: # FIX ME
            self.last_reaction = self.last_event
        else:
            self.last_reaction = None

        if __debug__:
#             log.info('\n%d: t=%g dt=%g\nevent=%s reactions=%d rejectedmoves=%d' % (self.stepCount, self.t, self.dt, self.last_event, self.reactionEvents, self.rejectedMoves))
            pass

        
        # Execute event.
        event.func(event.rr)

        if self.t != utils.INF: # and self.scheduler.getTopTime() == utils.INF
            nextTime = self.get_next_time()
            self.dt = nextTime - self.t
        else:
            self.dt = 0.0 # inf - inf == nan

    def fire(self, rr):
        for id in rr.reactants:
            st = self.model.get_species_type_by_id(id)
            self.remove_particles(st, 1)

        for id in rr.products:
            st = self.model.get_species_type_by_id(id)
            self.throw_in_particles(st, 1)

        for rr2 in self.dependencies[rr]:
            event = self.scheduler[rr2.eventID]
            dt = self.get_step_interval(rr2)
            self.update_event_time(self.t + dt, event)

        self.add_reaction_event(rr)

    def update(self, rr):
        self.dependencies.pop(rr)
        self.add_reaction_event(rr)

#     def remove_particles(self, st, n):
#         GillespieSimulatorBase.remove_particles(self, st.id, n)

    def throw_in_particles(self, st, n, surface=None):
        if not self.speciesDict.has_key(st.id):
            self.add_species(st.id)

        GillespieSimulatorBase.throw_in_particles(self, st, n)

    def get_species_id(self):
        return self.speciesDict.keys()

    def get_species(self):
        return []

    def add_species(self, id):
        GillespieSimulatorBase.add_species(self, id)

        st1 = self.model.get_species_type_by_id(id)
        rules = self.get_reaction_rule1(st1)
        for rr in rules:
            self.add_update_event(rr)

        for id2 in self.speciesDict.keys():
            st2 = self.model.get_species_type_by_id(id2)
            rules = self.get_reaction_rule2(st1, st2)
            for rr in rules:
                self.add_update_event(rr)
        
    def add_event(self, t, func, arg):
        return self.scheduler.add(GillespieEvent(t, func, arg))

    def add_reaction_event(self, rr):
        rr.eventType = EventType.SINGLE_REACTION # FIX ME

        dt = self.get_step_interval(rr)
        eventID = self.add_event(self.t + dt,
                                 Delegate(self, GillespieSimulator.fire),
                                 rr)
        if __debug__:
            log.info('addReactionEvent: #%d (t=%g)' % (eventID, self.t + dt))

        rr.eventID = eventID

        self.update_event_dependency(rr)

        for id, event in self.scheduler:
            rr2 = event.rr
            if rr == rr2:
                continue

            if rr.is_dependent_on(rr2):
                self.dependencies[rr2].append(rr)

    def add_update_event(self, rr):
        rr.eventType = EventType.SINGLE_ESCAPE # FIX ME

        eventID = self.add_event(self.t,
                                 Delegate(self, GillespieSimulator.update),
                                 rr)
        if __debug__:
            log.info('addUpdateEvent: #%d (t=%g)' % (eventID, self.t))

        rr.eventID = eventID
        self.dependencies[rr] = []

    def remove_event(self, rr):
        if __debug__:
            log.info('removeEvent: #%d' % rr.eventID)

        self.scheduler.removeEvent(rr.eventID)
        rr.eventID = None
        rr.eventType = None

    def update_event_time(self, t, event):
        if __debug__:
            log.info('updateEventTime: #%d (t=%g)' % (event.rr.eventID, t))

        self.scheduler.update((event.rr.eventID,
                               GillespieEvent(t, event.func, event.rr)))

#     def update_all_event_time(self):
#         for i in range(self.scheduler.size):
#             event = self.scheduler.getEventByIndex(i)
#             rr = event.getArg()

#             if rr.eventType == EventType.SINGLE_REACTION: # FIX ME
#                 dt = self.get_step_interval(rr)
#                 self.update_event_time(self.t + dt, rr)
#             else:
#                 assert self.t == event.getTime()
#                 self.update_event_time(self.t, rr)

    def update_event_dependency(self, rr1):
        self.dependencies[rr1] = []

        for id, event in self.scheduler:
            rr2 = event.rr
            if rr1 == rr2:
                continue

            if rr2.is_dependent_on(rr1):
                self.dependencies[rr1].append(rr2)

        self.dependencies[rr1].sort()

    def update_all_event_dependency(self):
        self.dependencies = {}

        for i in range(self.scheduler.size):
            rr = self.scheduler.getEventByIndex(i).getArg()
            self.update_event_dependency(rr)

    def interrupted(self, rr):
        '''return bool for efficiency.
        '''
        if float(rr.k) == 0.0:
            return False

        interrupt = False

        for id in rr.reactants:
            if self.govern(id):
                st = self.model.get_species_type_by_id(id)
                self.remove_particles(st, 1)
                interrupt = True

        for id in rr.products:
            if self.govern(id):
                st = self.model.get_species_type_by_id(id)
                self.throw_in_particles(st, 1)
                interrupt = True

        if not interrupt:
            return False

        rr1 = self.create_reaction_rule_cache(rr)
        for i in range(self.scheduler.size):
            event = self.scheduler.getEventByIndex(i)
            rr2 = event.getArg()
            if rr2.is_dependent_on(rr1):
                dt = self.get_step_interval(rr2)
                self.update_event_time(self.t + dt, event)

        dt = self.get_next_time() - self.t
        self.dt = dt

        return True

    def check(self):
        GillespieSimulatorBase.check()

        assert self.scheduler.check()
        assert self.t >= 0.0

    def dump_scheduler(self):
        for i in range(self.scheduler.size):
            event = self.scheduler.getEventByIndex(i)
            print i, event.getTime(), event.getArg()

    def dump(self):
        self.dump_scheduler()


if __name__ == '__main__':

    def main():
        pass


    main()

