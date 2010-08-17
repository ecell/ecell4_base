e754002b (shafi                   2008-04-27 06:07:06 +0000    1) #!/usr/env python
fb26ad8a (shafi                   2006-12-23 04:36:30 +0000    2) 
794888ab (shafi                   2008-04-29 01:31:15 +0000    3) 
9a603070 (shafi                   2009-11-26 06:11:41 +0000    4) from weakref import ref
fb26ad8a (shafi                   2006-12-23 04:36:30 +0000    5) import math
fb26ad8a (shafi                   2006-12-23 04:36:30 +0000    6) 
fb26ad8a (shafi                   2006-12-23 04:36:30 +0000    7) import numpy
fb26ad8a (shafi                   2006-12-23 04:36:30 +0000    8) 
109ada6f (Moriyoshi Koizumi       2010-02-02 16:44:13 +0900    9) from _gfrd import (
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900   10)     Event,
109ada6f (Moriyoshi Koizumi       2010-02-02 16:44:13 +0900   11)     EventScheduler,
109ada6f (Moriyoshi Koizumi       2010-02-02 16:44:13 +0900   12)     Particle,
109ada6f (Moriyoshi Koizumi       2010-02-02 16:44:13 +0900   13)     SphericalShell,
109ada6f (Moriyoshi Koizumi       2010-02-02 16:44:13 +0900   14)     SphericalShellContainer,
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900   15)     CylindricalShellContainer,
109ada6f (Moriyoshi Koizumi       2010-02-02 16:44:13 +0900   16)     DomainIDGenerator,
109ada6f (Moriyoshi Koizumi       2010-02-02 16:44:13 +0900   17)     ShellIDGenerator,
109ada6f (Moriyoshi Koizumi       2010-02-02 16:44:13 +0900   18)     DomainID,
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   19)     ParticleContainer,
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   20)     CuboidalRegion,
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   21)     CylindricalSurface,
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   22)     PlanarSurface,
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   23)     _random_vector
109ada6f (Moriyoshi Koizumi       2010-02-02 16:44:13 +0900   24)     )
fb26ad8a (shafi                   2006-12-23 04:36:30 +0000   25) 
f2a9c339 (thomie                  2010-05-12 03:48:34 +0200   26) from _greens_functions import EventType
f2a9c339 (thomie                  2010-05-12 03:48:34 +0200   27) 
a7e9eaf5 (shafi                   2007-01-12 21:49:06 +0000   28) from gfrdbase import *
72d01cc4 (thomie                  2010-02-03 17:33:23 +0900   29) from single import *
72d01cc4 (thomie                  2010-02-03 17:33:23 +0900   30) from pair import *
72d01cc4 (thomie                  2010-02-03 17:33:23 +0900   31) from multi import *
3f193d59 (moriyoshi               2009-11-27 07:38:17 +0000   32) from utils import *
fb26ad8a (shafi                   2006-12-23 04:36:30 +0000   33) 
44f4f9b4 (moriyoshi               2009-04-27 08:02:34 +0000   34) import logging
5069c176 (moriyoshi               2009-12-09 05:53:26 +0000   35) import os
44f4f9b4 (moriyoshi               2009-04-27 08:02:34 +0000   36) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900   37) from bd import DEFAULT_DT_FACTOR
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900   38) 
596a0d04 (shafi                   2009-12-18 08:03:33 +0000   39) log = logging.getLogger('ecell')
794888ab (shafi                   2008-04-29 01:31:15 +0000   40) 
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   41) def create_default_single(domain_id, pid_particle_pair, shell_id, rt, surface):
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   42)     if isinstance(surface, CuboidalRegion):
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   43)         return SphericalSingle(domain_id, pid_particle_pair, shell_id, rt, surface)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   44)     elif isinstance(surface, CylindricalSurface):
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   45)         return CylindricalSurfaceSingle(domain_id, pid_particle_pair, shell_id, rt, surface)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   46)     elif isinstance(surface, PlanarSurface):
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   47)         return PlanarSurfaceSingle(domain_id, pid_particle_pair, shell_id, rt, surface)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   48) 
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   49) def create_default_pair(domain_id, com, single1, single2, shell_id, 
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   50)                         r0, shell_size, rt, surface):
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   51)     if isinstance(surface, CuboidalRegion):
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   52)         return SphericalPair(domain_id, com, single1, single2, shell_id, r0, shell_size, rt, surface)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   53)     elif isinstance(surface, CylindricalSurface):
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   54)         return CylindricalSurfacePair(domain_id, com, single1, single2, shell_id, r0, shell_size, rt, surface)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   55)     elif isinstance(surface, PlanarSurface):
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   56)         return PlanarSurfacePair(domain_id, com, single1, single2, shell_id, r0, shell_size, rt, surface)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   57) 
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   58) 
986c2370 (Moriyoshi Koizumi       2010-07-23 20:03:14 +0900   59) class DomainEvent(Event):
986c2370 (Moriyoshi Koizumi       2010-07-23 20:03:14 +0900   60)     __slot__ = ['data']
986c2370 (Moriyoshi Koizumi       2010-07-23 20:03:14 +0900   61)     def __init__(self, time, domain):
986c2370 (Moriyoshi Koizumi       2010-07-23 20:03:14 +0900   62)         Event.__init__(self, time)
986c2370 (Moriyoshi Koizumi       2010-07-23 20:03:14 +0900   63)         self.data = domain
986c2370 (Moriyoshi Koizumi       2010-07-23 20:03:14 +0900   64) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900   65) class Delegate(object):
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900   66)     def __init__(self, obj, method, arg):
a034b53b (thomie                  2010-03-19 20:09:23 +0900   67)         self.ref = ref(obj)
1b002b58 (shafi                   2007-10-31 13:19:07 +0000   68)         self.method = method
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900   69)         self.arg = arg
1b002b58 (shafi                   2007-10-31 13:19:07 +0000   70) 
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900   71)     def __call__(self):
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900   72)         return self.method(self.ref(), self.arg)
1b002b58 (shafi                   2007-10-31 13:19:07 +0000   73) 
9a603070 (shafi                   2009-11-26 06:11:41 +0000   74) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900   75) class EGFRDSimulator(ParticleSimulatorBase):
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   76)     def __init__(self, world, rng, network_rules):
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900   77)         ParticleSimulatorBase.__init__(self, world, rng, network_rules)
acd00e2c (Moriyoshi Koizumi       2010-03-11 14:21:35 +0900   78) 
7790978e (thomie                  2010-03-19 20:27:36 +0900   79)         self.domain_id_generator = DomainIDGenerator(0)
7790978e (thomie                  2010-03-19 20:27:36 +0900   80)         self.shell_id_generator = ShellIDGenerator(0)
f87104f5 (shafi                   2008-06-08 04:12:33 +0000   81) 
2e23fcbb (shafi                   2008-07-17 22:12:24 +0000   82)         self.MULTI_SHELL_FACTOR = 0.05
c0dffe52 (shafi                   2008-07-18 19:39:47 +0000   83)         self.SINGLE_SHELL_FACTOR = 0.1
26003e3b (shafi                   2008-04-22 08:43:35 +0000   84) 
7790978e (thomie                  2010-03-19 20:27:36 +0900   85)         self.is_dirty = True
4c23d166 (shafi                   2007-03-28 07:17:41 +0000   86)         self.scheduler = EventScheduler()
09b80656 (shafi                   2007-02-06 02:30:46 +0000   87) 
7790978e (thomie                  2010-03-19 20:27:36 +0900   88)         self.user_max_shell_size = numpy.inf
a7e32062 (shafi                   2007-09-08 02:27:45 +0000   89) 
2c5614ad (shafi                   2009-12-15 04:46:20 +0000   90)         self.domains = {}
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000   91) 
edfb1a35 (shafi                   2007-12-12 02:15:51 +0000   92)         self.reset()
a7e32062 (shafi                   2007-09-08 02:27:45 +0000   93) 
7790978e (thomie                  2010-03-19 20:27:36 +0900   94)     def get_matrix_cell_size(self):
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900   95)         return self.containers[0].cell_size
b99a9cf5 (shafi                   2007-10-31 11:45:54 +0000   96) 
7790978e (thomie                  2010-03-19 20:27:36 +0900   97)     def get_next_time(self):
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900   98)         if self.scheduler.size == 0:
0dca5c43 (shafi                   2008-07-17 18:01:08 +0000   99)             return self.t
0dca5c43 (shafi                   2008-07-17 18:01:08 +0000  100) 
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  101)         return self.scheduler.top[1].time
2c6d4abb (shafi                   2008-03-14 21:01:21 +0000  102) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  103)     def set_user_max_shell_size(self, size):
7790978e (thomie                  2010-03-19 20:27:36 +0900  104)         self.user_max_shell_size = size
11c410e8 (shafi                   2008-06-09 01:10:56 +0000  105) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  106)     def get_user_max_shell_size(self):
7790978e (thomie                  2010-03-19 20:27:36 +0900  107)         return self.user_max_shell_size
2ced4967 (shafi                   2007-10-17 11:20:20 +0000  108) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  109)     def get_max_shell_size(self):
7790978e (thomie                  2010-03-19 20:27:36 +0900  110)         return min(self.get_matrix_cell_size() * .5 / SAFETY,
7790978e (thomie                  2010-03-19 20:27:36 +0900  111)                    self.user_max_shell_size)
2ced4967 (shafi                   2007-10-17 11:20:20 +0000  112) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900  113)     def reset(self):
edfb1a35 (shafi                   2007-12-12 02:15:51 +0000  114)         self.t = 0.0
89d8bae4 (shafi                   2008-06-23 05:53:56 +0000  115)         self.dt = 0.0
7790978e (thomie                  2010-03-19 20:27:36 +0900  116)         self.step_counter = 0
616a0341 (thomie                  2010-02-22 17:13:13 +0900  117)         self.single_steps = {EventType.SINGLE_ESCAPE:0,
616a0341 (thomie                  2010-02-22 17:13:13 +0900  118)                              EventType.SINGLE_REACTION:0}
616a0341 (thomie                  2010-02-22 17:13:13 +0900  119)         self.pair_steps = {EventType.SINGLE_REACTION:0,
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  120)                            EventType.IV_REACTION:0,
616a0341 (thomie                  2010-02-22 17:13:13 +0900  121)                            EventType.IV_ESCAPE:0,
616a0341 (thomie                  2010-02-22 17:13:13 +0900  122)                            EventType.COM_ESCAPE:0}
616a0341 (thomie                  2010-02-22 17:13:13 +0900  123)         self.multi_steps = {EventType.MULTI_ESCAPE:0,
616a0341 (thomie                  2010-02-22 17:13:13 +0900  124)                             EventType.MULTI_REACTION:0, 2:0}
7790978e (thomie                  2010-03-19 20:27:36 +0900  125)         self.zero_steps = 0
7790978e (thomie                  2010-03-19 20:27:36 +0900  126)         self.rejected_moves = 0
7790978e (thomie                  2010-03-19 20:27:36 +0900  127)         self.reaction_events = 0
7790978e (thomie                  2010-03-19 20:27:36 +0900  128)         self.last_event = None
7790978e (thomie                  2010-03-19 20:27:36 +0900  129)         self.last_reaction = None
d99383b0 (shafi                   2008-01-10 05:50:40 +0000  130) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  131)         self.is_dirty = True
edfb1a35 (shafi                   2007-12-12 02:15:51 +0000  132) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900  133)     def initialize(self):
a034b53b (thomie                  2010-03-19 20:09:23 +0900  134)         ParticleSimulatorBase.initialize(self)
e46219c4 (shafi                   2008-03-21 00:31:10 +0000  135) 
b2c8b16f (shafi                   2007-02-10 01:11:11 +0000  136)         self.scheduler.clear()
acd00e2c (Moriyoshi Koizumi       2010-03-11 14:21:35 +0900  137)         self.containers = [SphericalShellContainer(self.world.world_size, 
acd00e2c (Moriyoshi Koizumi       2010-03-11 14:21:35 +0900  138)                                                    self.world.matrix_size),
acd00e2c (Moriyoshi Koizumi       2010-03-11 14:21:35 +0900  139)                            CylindricalShellContainer(self.world.world_size, 
acd00e2c (Moriyoshi Koizumi       2010-03-11 14:21:35 +0900  140)                                                      self.world.matrix_size)]
61d38d3d (thomie                  2010-02-04 21:29:13 +0900  141)         self.domains = {}
b2c8b16f (shafi                   2007-02-10 01:11:11 +0000  142) 
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  143)         singles = []
3ef80c8f (thomie                  2010-08-09 11:04:51 +0200  144) 
3ef80c8f (thomie                  2010-08-09 11:04:51 +0200  145)         # Fix order of adding particles (always, or at least in debug mode).
3ef80c8f (thomie                  2010-08-09 11:04:51 +0200  146)         pid_particle_pairs = list(self.world)
3ef80c8f (thomie                  2010-08-09 11:04:51 +0200  147)         pid_particle_pairs.sort()
3ef80c8f (thomie                  2010-08-09 11:04:51 +0200  148) 
3ef80c8f (thomie                  2010-08-09 11:04:51 +0200  149)         for pid_particle_pair in pid_particle_pairs:
7790978e (thomie                  2010-03-19 20:27:36 +0900  150)             single = self.create_single(pid_particle_pair)
69126931 (Moriyoshi Koizumi       2010-03-15 16:56:22 +0900  151)             if __debug__:
69126931 (Moriyoshi Koizumi       2010-03-15 16:56:22 +0900  152)                 log.debug("%s as single %s", pid_particle_pair[0], single.domain_id)
69126931 (Moriyoshi Koizumi       2010-03-15 16:56:22 +0900  153)             singles.append(single)
acd00e2c (Moriyoshi Koizumi       2010-03-11 14:21:35 +0900  154)         assert len(singles) == self.world.num_particles
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  155)         for single in singles:
7790978e (thomie                  2010-03-19 20:27:36 +0900  156)             self.add_single_event(single)
5bf2f937 (shafi                   2007-02-08 03:34:33 +0000  157) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  158)         self.is_dirty = False
09b80656 (shafi                   2007-02-06 02:30:46 +0000  159) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900  160)     def stop(self, t):
049b6bdd (mozo                    2009-01-26 06:54:38 +0000  161)         if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900  162)             log.info('stop at %g' % t)
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  163) 
d1b2fe72 (shafi                   2007-10-30 08:45:30 +0000  164)         if self.t == t:
d1b2fe72 (shafi                   2007-10-30 08:45:30 +0000  165)             return
d1b2fe72 (shafi                   2007-10-30 08:45:30 +0000  166) 
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  167)         if t >= self.scheduler.top[1].time:
68760ef5 (shafi                   2008-06-22 02:16:31 +0000  168)             raise RuntimeError, 'Stop time >= next event time.'
b8ed0d73 (shafi                   2007-10-25 09:34:49 +0000  169) 
725f9edd (shafi                   2008-07-19 00:39:16 +0000  170)         if t < self.t:
725f9edd (shafi                   2008-07-19 00:39:16 +0000  171)             raise RuntimeError, 'Stop time < current time.'
725f9edd (shafi                   2008-07-19 00:39:16 +0000  172) 
fec1b911 (shafi                   2007-07-09 07:26:59 +0000  173)         self.t = t
68c4d299 (shafi@tcs1.gsc.riken.jp 2010-02-03 13:24:14 +0900  174) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  175)         non_single_list = []
fec1b911 (shafi                   2007-07-09 07:26:59 +0000  176) 
fec1b911 (shafi                   2007-07-09 07:26:59 +0000  177)         # first burst all Singles.
bed38b09 (Moriyoshi Koizumi       2010-07-05 16:21:47 +0900  178)         for id, event in self.scheduler:
bed38b09 (Moriyoshi Koizumi       2010-07-05 16:21:47 +0900  179)             obj = event.data
a034b53b (thomie                  2010-03-19 20:09:23 +0900  180)             if isinstance(obj, Pair) or isinstance(obj, Multi):
7790978e (thomie                  2010-03-19 20:27:36 +0900  181)                 non_single_list.append(obj)
a034b53b (thomie                  2010-03-19 20:09:23 +0900  182)             elif isinstance(obj, Single):
049b6bdd (mozo                    2009-01-26 06:54:38 +0000  183)                 if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900  184)                     log.debug('burst %s, last_time= %g' % 
7790978e (thomie                  2010-03-19 20:27:36 +0900  185)                           (str(obj), obj.last_time))
7790978e (thomie                  2010-03-19 20:27:36 +0900  186)                 self.burst_single(obj)
5ba6f664 (shafi                   2008-03-26 23:47:39 +0000  187)             else:
faba9096 (shafi                   2008-04-22 05:17:09 +0000  188)                 assert False, 'do not reach here'
5ba6f664 (shafi                   2008-03-26 23:47:39 +0000  189) 
fec1b911 (shafi                   2007-07-09 07:26:59 +0000  190) 
8b0df249 (shafi                   2008-06-20 05:28:58 +0000  191)         # then burst all Pairs and Multis.
049b6bdd (mozo                    2009-01-26 06:54:38 +0000  192)         if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900  193)             log.debug('burst %s' % non_single_list)
7790978e (thomie                  2010-03-19 20:27:36 +0900  194)         self.burst_objs(non_single_list)
fec1b911 (shafi                   2007-07-09 07:26:59 +0000  195) 
fec1b911 (shafi                   2007-07-09 07:26:59 +0000  196)         self.dt = 0.0
fec1b911 (shafi                   2007-07-09 07:26:59 +0000  197) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900  198)     def step(self):
7790978e (thomie                  2010-03-19 20:27:36 +0900  199)         self.last_reaction = None
8cbd82e3 (shafi                   2007-07-10 07:16:21 +0000  200) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  201)         if self.is_dirty:
5bf2f937 (shafi                   2007-02-08 03:34:33 +0000  202)             self.initialize()
36086d2b (shafi                   2008-07-03 10:54:30 +0000  203)             
5069c176 (moriyoshi               2009-12-09 05:53:26 +0000  204)         if __debug__:
fd25dc19 (thomie                  2010-03-12 13:03:35 +0900  205)             if int("0" + os.environ.get("ECELL_CHECK", ""), 10):
fd25dc19 (thomie                  2010-03-12 13:03:35 +0900  206)                 self.check()
f87104f5 (shafi                   2008-06-08 04:12:33 +0000  207)         
7790978e (thomie                  2010-03-19 20:27:36 +0900  208)         self.step_counter += 1
5a5dacbb (shafi                   2007-09-24 23:14:31 +0000  209) 
e463405f (thomie                  2010-03-01 10:54:33 +0900  210)         if __debug__:
34bb0cae (thomie                  2010-06-30 15:55:51 +0200  211)             if self.scheduler.size == 0:
e463405f (thomie                  2010-03-01 10:54:33 +0900  212)                 raise RuntimeError('No particles in scheduler.')
e463405f (thomie                  2010-03-01 10:54:33 +0900  213) 
3ef80c8f (thomie                  2010-08-09 11:04:51 +0200  214)         id, event = self.scheduler.pop()
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  215)         self.t, self.last_event = event.time, event
41aaa33c (shafi                   2007-02-12 01:35:08 +0000  216) 
049b6bdd (mozo                    2009-01-26 06:54:38 +0000  217)         if __debug__:
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000  218)             domain_counts = self.count_domains()
a034b53b (thomie                  2010-03-19 20:09:23 +0900  219)             log.info('\n%d: t=%g dt=%g\tSingles: %d, Pairs: %d, Multis: %d'
7790978e (thomie                  2010-03-19 20:27:36 +0900  220)                      % ((self.step_counter, self.t, self.dt) + domain_counts))
3ef80c8f (thomie                  2010-08-09 11:04:51 +0200  221)             log.info('event=%d reactions=%d rejectedmoves=%d' 
3ef80c8f (thomie                  2010-08-09 11:04:51 +0200  222)                      % (id, self.reaction_events, 
7790978e (thomie                  2010-03-19 20:27:36 +0900  223)                         self.rejected_moves))
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900  224)        
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900  225)         for klass, f in self.dispatch:
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900  226)             if isinstance(event.data, klass):
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900  227)                 f(self, event.data)
41aaa33c (shafi                   2007-02-12 01:35:08 +0000  228) 
143ba3f0 (thomie                  2010-02-23 15:44:44 +0900  229)         if __debug__:
34bb0cae (thomie                  2010-06-30 15:55:51 +0200  230)             if self.scheduler.size == 0:
143ba3f0 (thomie                  2010-02-23 15:44:44 +0900  231)                 raise RuntimeError('Zero particles left.')
143ba3f0 (thomie                  2010-02-23 15:44:44 +0900  232) 
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  233)         next_time = self.scheduler.top[1].time
7790978e (thomie                  2010-03-19 20:27:36 +0900  234)         self.dt = next_time - self.t
df639b72 (shafi                   2007-06-29 04:00:53 +0000  235) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000  236)         # assert if not too many successive dt=0 steps occur.
d29ea599 (shafi                   2009-04-14 05:18:15 +0000  237)         if __debug__:
d29ea599 (shafi                   2009-04-14 05:18:15 +0000  238)             if self.dt == 0:
7790978e (thomie                  2010-03-19 20:27:36 +0900  239)                 self.zero_steps += 1
34bb0cae (thomie                  2010-06-30 15:55:51 +0200  240)                 if self.zero_steps >= max(self.scheduler.size * 3, 10):
d29ea599 (shafi                   2009-04-14 05:18:15 +0000  241)                     raise RuntimeError, 'too many dt=zero steps.  simulator halted?'
d29ea599 (shafi                   2009-04-14 05:18:15 +0000  242)             else:
7790978e (thomie                  2010-03-19 20:27:36 +0900  243)                 self.zero_steps = 0
8b2a7032 (shafi                   2008-02-15 04:02:20 +0000  244) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  245)     def create_single(self, pid_particle_pair):
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  246)         rt = self.network_rules.query_reaction_rule(pid_particle_pair[1].sid)
7790978e (thomie                  2010-03-19 20:27:36 +0900  247)         domain_id = self.domain_id_generator()
7790978e (thomie                  2010-03-19 20:27:36 +0900  248)         shell_id = self.shell_id_generator()
df083205 (shafi                   2007-09-19 23:30:15 +0000  249) 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  250)         # Get surface.
acd00e2c (Moriyoshi Koizumi       2010-03-11 14:21:35 +0900  251)         species = self.world.get_species(pid_particle_pair[1].sid)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  252)         surface = self.world.get_structure(species.structure_id)
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  253) 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  254)         # Create single. The type of the single that will be created depends 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  255)         # on the surface this particle is on. Either SphericalSingle, 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  256)         # PlanarSurfaceSingle, or CylindricalSurfaceSingle.
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  257)         single = create_default_single(domain_id, pid_particle_pair, shell_id, rt, surface)
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  258) 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  259)         single.initialize(self.t)
7790978e (thomie                  2010-03-19 20:27:36 +0900  260)         self.move_shell(single.shell_id_shell_pair)
3986c94d (thomie                  2010-02-05 22:29:38 +0900  261)         self.domains[domain_id] = single
3986c94d (thomie                  2010-02-05 22:29:38 +0900  262)         return single
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  263) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  264)     def create_pair(self, single1, single2, com, r0, shell_size):
df083205 (shafi                   2007-09-19 23:30:15 +0000  265)         assert single1.dt == 0.0
df083205 (shafi                   2007-09-19 23:30:15 +0000  266)         assert single2.dt == 0.0
df083205 (shafi                   2007-09-19 23:30:15 +0000  267) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  268)         rt = self.network_rules.query_reaction_rule(single1.pid_particle_pair[1].sid, single2.pid_particle_pair[1].sid)[0]
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000  269) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  270)         domain_id = self.domain_id_generator()
7790978e (thomie                  2010-03-19 20:27:36 +0900  271)         shell_id = self.shell_id_generator()
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  272) 
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  273)         pos1 = single1.shell.shape.position
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  274)         pos2 = single2.shell.shape.position
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  275) 
a60a916f (thomie                  2010-02-22 13:55:57 +0900  276)         # Get surface.
acd00e2c (Moriyoshi Koizumi       2010-03-11 14:21:35 +0900  277)         species = self.world.get_species(single1.pid_particle_pair[1].sid)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  278)         surface = self.world.get_structure(species.structure_id)
033e5a59 (shafi                   2009-12-21 07:12:54 +0000  279) 
a60a916f (thomie                  2010-02-22 13:55:57 +0900  280)         # Create pair. The type of the pair that will be created depends on 
a60a916f (thomie                  2010-02-22 13:55:57 +0900  281)         # the surface the particles are on. Either SphericalPair, 
a60a916f (thomie                  2010-02-22 13:55:57 +0900  282)         # PlanarSurfacePair, or CylindricalSurfacePair.
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  283)         pair = create_default_pair(domain_id, com, single1, single2, shell_id, 
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  284)                                    r0, shell_size, rt, surface)
033e5a59 (shafi                   2009-12-21 07:12:54 +0000  285) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900  286)         pair.initialize(self.t)
033e5a59 (shafi                   2009-12-21 07:12:54 +0000  287) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  288)         self.move_shell(pair.shell_id_shell_pair)
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000  289)         self.domains[domain_id] = pair
df083205 (shafi                   2007-09-19 23:30:15 +0000  290)         return pair
033e5a59 (shafi                   2009-12-21 07:12:54 +0000  291) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  292)     def create_multi(self):
7790978e (thomie                  2010-03-19 20:27:36 +0900  293)         domain_id = self.domain_id_generator()
d2feaaec (thomie                  2010-03-01 12:28:04 +0900  294)         if __debug__:
d2feaaec (thomie                  2010-03-01 12:28:04 +0900  295)             try:
d2feaaec (thomie                  2010-03-01 12:28:04 +0900  296)                 # Option to make multis run faster for nicer visualization.
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  297)                 dt_factor = DEFAULT_DT_FACTOR * self.bd_dt_factor
d2feaaec (thomie                  2010-03-01 12:28:04 +0900  298)             except AttributeError:
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  299)                 dt_factor = DEFAULT_DT_FACTOR 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  300)         else:
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  301)             dt_factor = DEFAULT_DT_FACTOR
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  302)         multi = Multi(domain_id, self, dt_factor)
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  303)         self.domains[domain_id] = multi
faba9096 (shafi                   2008-04-22 05:17:09 +0000  304)         return multi
faba9096 (shafi                   2008-04-22 05:17:09 +0000  305) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  306)     def move_single(self, single, position, radius=None):
7790978e (thomie                  2010-03-19 20:27:36 +0900  307)         self.move_single_shell(single, position, radius)
7790978e (thomie                  2010-03-19 20:27:36 +0900  308)         self.move_single_particle(single, position)
3986c94d (thomie                  2010-02-05 22:29:38 +0900  309) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  310)     def move_single_shell(self, single, position, radius=None):
fb5d7555 (thomie                  2010-02-26 18:42:41 +0900  311)         if radius == None:
fb5d7555 (thomie                  2010-02-26 18:42:41 +0900  312)             # By default, don't change radius.
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  313)             radius = single.shell.shape.radius
3986c94d (thomie                  2010-02-05 22:29:38 +0900  314) 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  315)         # Reuse shell_id and domain_id.
e9d894b2 (thomie                  2010-02-26 21:42:07 +0900  316)         shell_id = single.shell_id
3986c94d (thomie                  2010-02-05 22:29:38 +0900  317)         domain_id = single.domain_id
3986c94d (thomie                  2010-02-05 22:29:38 +0900  318) 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  319)         # Replace shell.
7790978e (thomie                  2010-03-19 20:27:36 +0900  320)         shell = single.create_new_shell(position, radius, domain_id)
e9d894b2 (thomie                  2010-02-26 21:42:07 +0900  321)         shell_id_shell_pair = (shell_id, shell) 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  322) 
e9d894b2 (thomie                  2010-02-26 21:42:07 +0900  323)         single.shell_id_shell_pair = shell_id_shell_pair
7790978e (thomie                  2010-03-19 20:27:36 +0900  324)         self.move_shell(shell_id_shell_pair)
3986c94d (thomie                  2010-02-05 22:29:38 +0900  325) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  326)     def move_single_particle(self, single, position):
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000  327)         new_pid_particle_pair = (single.pid_particle_pair[0],
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000  328)                           Particle(position,
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000  329)                                    single.pid_particle_pair[1].radius,
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000  330)                                    single.pid_particle_pair[1].D,
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000  331)                                    single.pid_particle_pair[1].sid))
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000  332)         single.pid_particle_pair = new_pid_particle_pair
3986c94d (thomie                  2010-02-05 22:29:38 +0900  333) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  334)         self.world.update_particle(new_pid_particle_pair)
8de778cb (shafi                   2008-07-22 00:28:00 +0000  335) 
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900  336)     def get_container(self, shell):
737c6f80 (thomie                  2010-03-13 00:19:25 +0900  337)         if type(shell) is SphericalShell:
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900  338)             return self.containers[0]
737c6f80 (thomie                  2010-03-13 00:19:25 +0900  339)         elif type(shell) is CylindricalShell:
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900  340)             return self.containers[1]
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900  341) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  342)     def remove_domain(self, obj):
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  343)         if __debug__:
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  344)             log.info("remove_domain: %s" % obj)
74304fc9 (moriyoshi               2009-12-09 07:15:51 +0000  345)         del self.domains[obj.domain_id]
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900  346)         for shell_id, shell in obj.shell_list:
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900  347)             container = self.get_container(shell)
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900  348)             del container[shell_id]
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900  349) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  350)     def move_shell(self, shell_id_shell_pair):
e9d894b2 (thomie                  2010-02-26 21:42:07 +0900  351)         shell = shell_id_shell_pair[1]
e9d894b2 (thomie                  2010-02-26 21:42:07 +0900  352)         container = self.get_container(shell)
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900  353)         container.update(shell_id_shell_pair)
8de778cb (shafi                   2008-07-22 00:28:00 +0000  354) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  355)     def add_single_event(self, single):
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  356)         event_id = self.scheduler.add(
986c2370 (Moriyoshi Koizumi       2010-07-23 20:03:14 +0900  357)             DomainEvent(self.t + single.dt, single))
ebb84956 (moriyoshi               2009-04-15 07:44:25 +0000  358)         if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900  359)             log.info('add_single_event: #%d (t=%g)' % (
7790978e (thomie                  2010-03-19 20:27:36 +0900  360)                event_id, self.t + single.dt))
7790978e (thomie                  2010-03-19 20:27:36 +0900  361)         single.event_id = event_id
dd6ae27b (shafi                   2007-11-27 00:09:18 +0000  362) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  363)     def add_pair_event(self, pair):
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  364)         event_id = self.scheduler.add(
986c2370 (Moriyoshi Koizumi       2010-07-23 20:03:14 +0900  365)             DomainEvent(self.t + pair.dt, pair))
ebb84956 (moriyoshi               2009-04-15 07:44:25 +0000  366)         if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900  367)             log.info('add_pair_event: #%d (t=%g)' % (
7790978e (thomie                  2010-03-19 20:27:36 +0900  368)                event_id, self.t + pair.dt))
7790978e (thomie                  2010-03-19 20:27:36 +0900  369)         pair.event_id = event_id
b5dbc44d (shafi                   2007-04-13 04:21:54 +0000  370) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  371)     def add_multi_event(self, multi):
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  372)         event_id = self.scheduler.add(
986c2370 (Moriyoshi Koizumi       2010-07-23 20:03:14 +0900  373)             DomainEvent(self.t + multi.dt, multi))
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  374) 
ebb84956 (moriyoshi               2009-04-15 07:44:25 +0000  375)         if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900  376)             log.info('add_multi_event: #%d (t=%g)' % (
7790978e (thomie                  2010-03-19 20:27:36 +0900  377)                event_id, self.t + multi.dt))
7790978e (thomie                  2010-03-19 20:27:36 +0900  378)         multi.event_id = event_id
faba9096 (shafi                   2008-04-22 05:17:09 +0000  379) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900  380)     def removeEvent(self, event):
ebb84956 (moriyoshi               2009-04-15 07:44:25 +0000  381)         if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900  382)             log.info('removeEvent: #%d' % event.event_id)
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  383)         del self.scheduler[event.event_id]
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  384) 
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  385)     def update_single_event(self, t, single):
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  386)         if __debug__:
34bb0cae (thomie                  2010-06-30 15:55:51 +0200  387)             log.info('update_event: #%d (t=%g)' % (single.event_id, t))
986c2370 (Moriyoshi Koizumi       2010-07-23 20:03:14 +0900  388)         self.scheduler.update((single.event_id, DomainEvent(t, single)))
ba9b6458 (shafi                   2007-04-20 03:18:26 +0000  389) 
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  390)     def update_multi_event(self, t, multi):
ebb84956 (moriyoshi               2009-04-15 07:44:25 +0000  391)         if __debug__:
34bb0cae (thomie                  2010-06-30 15:55:51 +0200  392)             log.info('update_event: #%d (t=%g)' % (multi.event_id, t))
986c2370 (Moriyoshi Koizumi       2010-07-23 20:03:14 +0900  393)         self.scheduler.update((multi.event_id, DomainEvent(t, multi)))
b5dbc44d (shafi                   2007-04-13 04:21:54 +0000  394) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  395)     def burst_obj(self, obj):
049b6bdd (mozo                    2009-01-26 06:54:38 +0000  396)         if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900  397)             log.info('burst_obj: bursting %s' % obj)
794888ab (shafi                   2008-04-29 01:31:15 +0000  398) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900  399)         if isinstance(obj, Single):
3986c94d (thomie                  2010-02-05 22:29:38 +0900  400)             # TODO. Compare with gfrd.
7790978e (thomie                  2010-03-19 20:27:36 +0900  401)             self.burst_single(obj)
a034b53b (thomie                  2010-03-19 20:09:23 +0900  402)             bursted = [obj, ]
a034b53b (thomie                  2010-03-19 20:09:23 +0900  403)         elif isinstance(obj, Pair):  # Pair
7790978e (thomie                  2010-03-19 20:27:36 +0900  404)             single1, single2 = self.burst_pair(obj)
7790978e (thomie                  2010-03-19 20:27:36 +0900  405)             # Don't schedule events in burst/propagate_pair, because scheduling 
7790978e (thomie                  2010-03-19 20:27:36 +0900  406)             # is different after a single reaction in fire_pair.
7790978e (thomie                  2010-03-19 20:27:36 +0900  407)             self.add_single_event(single1)
7790978e (thomie                  2010-03-19 20:27:36 +0900  408)             self.add_single_event(single2)
a034b53b (thomie                  2010-03-19 20:09:23 +0900  409)             self.removeEvent(obj)
a034b53b (thomie                  2010-03-19 20:09:23 +0900  410)             bursted = [single1, single2]
faba9096 (shafi                   2008-04-22 05:17:09 +0000  411)         else:  # Multi
7790978e (thomie                  2010-03-19 20:27:36 +0900  412)             bursted = self.burst_multi(obj)
a034b53b (thomie                  2010-03-19 20:09:23 +0900  413)             self.removeEvent(obj)
9d0566db (moriyoshi               2009-04-23 09:09:24 +0000  414) 
9d0566db (moriyoshi               2009-04-23 09:09:24 +0000  415)         if __debug__:
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  416)             log.info('burst_obj: bursted=%s' % ', '.join(str(i) for i in bursted))
9d0566db (moriyoshi               2009-04-23 09:09:24 +0000  417) 
9d0566db (moriyoshi               2009-04-23 09:09:24 +0000  418)         return bursted
faba9096 (shafi                   2008-04-22 05:17:09 +0000  419) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  420)     def burst_objs(self, objs):
faba9096 (shafi                   2008-04-22 05:17:09 +0000  421)         bursted = []
faba9096 (shafi                   2008-04-22 05:17:09 +0000  422)         for obj in objs:
7790978e (thomie                  2010-03-19 20:27:36 +0900  423)             b = self.burst_obj(obj)
a034b53b (thomie                  2010-03-19 20:09:23 +0900  424)             bursted.extend(b)
faba9096 (shafi                   2008-04-22 05:17:09 +0000  425) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000  426)         return bursted
faba9096 (shafi                   2008-04-22 05:17:09 +0000  427) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  428)     def clear_volume(self, pos, radius, ignore=[]):
7790978e (thomie                  2010-03-19 20:27:36 +0900  429)         neighbors = self.get_neighbors_within_radius_no_sort(pos, radius, ignore)
7790978e (thomie                  2010-03-19 20:27:36 +0900  430)         return self.burst_objs(neighbors)
faba9096 (shafi                   2008-04-22 05:17:09 +0000  431) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  432)     def burst_non_multis(self, neighbors):
faba9096 (shafi                   2008-04-22 05:17:09 +0000  433)         bursted = []
faba9096 (shafi                   2008-04-22 05:17:09 +0000  434) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000  435)         for obj in neighbors:
a034b53b (thomie                  2010-03-19 20:09:23 +0900  436)             if not isinstance(obj, Multi):
7790978e (thomie                  2010-03-19 20:27:36 +0900  437)                 b = self.burst_obj(obj)
a034b53b (thomie                  2010-03-19 20:09:23 +0900  438)                 bursted.extend(b)
faba9096 (shafi                   2008-04-22 05:17:09 +0000  439)             else:
a034b53b (thomie                  2010-03-19 20:09:23 +0900  440)                 bursted.append(obj)
faba9096 (shafi                   2008-04-22 05:17:09 +0000  441) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000  442)         return bursted
faba9096 (shafi                   2008-04-22 05:17:09 +0000  443) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  444)     def fire_single_reaction(self, single):
7790978e (thomie                  2010-03-19 20:27:36 +0900  445)         reactant_species_radius = single.pid_particle_pair[1].radius
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000  446)         oldpos = single.pid_particle_pair[1].position
7790978e (thomie                  2010-03-19 20:27:36 +0900  447)         current_surface = single.surface
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  448)         
7790978e (thomie                  2010-03-19 20:27:36 +0900  449)         rt = single.draw_reaction_rule()
8bfd4ae0 (shafi                   2007-11-29 09:29:45 +0000  450) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900  451)         if len(rt.products) == 0:
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  452)             
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  453)             self.world.remove_particle(single.pid_particle_pair[0])
73377946 (shafi                   2008-08-05 23:30:41 +0000  454) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  455)             self.last_reaction = (rt, (single.pid_particle_pair[1], None), [])
73377946 (shafi                   2008-08-05 23:30:41 +0000  456) 
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  457)             
a034b53b (thomie                  2010-03-19 20:09:23 +0900  458)         elif len(rt.products) == 1:
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  459)             
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  460)             product_species = self.world.get_species(rt.products[0])
f0458afe (shafi                   2007-09-26 02:45:44 +0000  461) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  462)             if reactant_species_radius < product_species.radius:
7790978e (thomie                  2010-03-19 20:27:36 +0900  463)                 self.clear_volume(oldpos, product_species.radius)
faba9096 (shafi                   2008-04-22 05:17:09 +0000  464) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  465)             if self.world.check_overlap((oldpos, product_species.radius),
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  466)                                         single.pid_particle_pair[0]):
049b6bdd (mozo                    2009-01-26 06:54:38 +0000  467)                 if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900  468)                     log.info('no space for product particle.')
4f1494cd (shafi                   2007-12-08 08:24:56 +0000  469)                 raise NoSpace()
36ed16fa (shafi                   2007-10-19 03:55:57 +0000  470) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  471)             self.world.remove_particle(single.pid_particle_pair[0])
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  472)             newparticle = self.world.new_particle(product_species.id, oldpos)
7790978e (thomie                  2010-03-19 20:27:36 +0900  473)             newsingle = self.create_single(newparticle)
7790978e (thomie                  2010-03-19 20:27:36 +0900  474)             self.add_single_event(newsingle)
73377946 (shafi                   2008-08-05 23:30:41 +0000  475) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  476)             self.last_reaction = (rt, (single.pid_particle_pair[1], None), [newparticle])
73377946 (shafi                   2008-08-05 23:30:41 +0000  477) 
049b6bdd (mozo                    2009-01-26 06:54:38 +0000  478)             if __debug__:
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  479)                 log.info('product = %s' % str(newsingle))
f0458afe (shafi                   2007-09-26 02:45:44 +0000  480) 
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  481)             
a034b53b (thomie                  2010-03-19 20:09:23 +0900  482)         elif len(rt.products) == 2:
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  483)             product_species1 = self.world.get_species(rt.products[0])
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  484)             product_species2 = self.world.get_species(rt.products[1])
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  485)             
7790978e (thomie                  2010-03-19 20:27:36 +0900  486)             D1 = product_species1.D
7790978e (thomie                  2010-03-19 20:27:36 +0900  487)             D2 = product_species2.D
f0458afe (shafi                   2007-09-26 02:45:44 +0000  488)             D12 = D1 + D2
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  489)             
7790978e (thomie                  2010-03-19 20:27:36 +0900  490)             particle_radius1 = product_species1.radius
7790978e (thomie                  2010-03-19 20:27:36 +0900  491)             particle_radius2 = product_species2.radius
7790978e (thomie                  2010-03-19 20:27:36 +0900  492)             particle_radius12 = particle_radius1 + particle_radius2
f0458afe (shafi                   2007-09-26 02:45:44 +0000  493) 
3c3e55d8 (shafi                   2008-02-19 00:44:24 +0000  494)             # clean up space.
7790978e (thomie                  2010-03-19 20:27:36 +0900  495)             rad = max(particle_radius12 * (D1 / D12) + particle_radius1,
7790978e (thomie                  2010-03-19 20:27:36 +0900  496)                       particle_radius12 * (D2 / D12) + particle_radius2)
3c3e55d8 (shafi                   2008-02-19 00:44:24 +0000  497) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  498)             self.clear_volume(oldpos, rad)
3c3e55d8 (shafi                   2008-02-19 00:44:24 +0000  499) 
525644ac (shafi@tcs1.gsc.riken.jp 2010-02-15 18:16:47 +0900  500)             for _ in range(self.dissociation_retry_moves):
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  501)                 vector = _random_vector(current_surface, particle_radius12 *
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  502)                                        MINIMAL_SEPERATION_FACTOR, self.rng)
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  503)             
f0458afe (shafi                   2007-09-26 02:45:44 +0000  504)                 # place particles according to the ratio D1:D2
f0458afe (shafi                   2007-09-26 02:45:44 +0000  505)                 # this way, species with D=0 doesn't move.
f0458afe (shafi                   2007-09-26 02:45:44 +0000  506)                 # FIXME: what if D1 == D2 == 0?
d522a419 (shafi                   2007-11-08 09:42:25 +0000  507) 
a41e0fdc (shafi                   2008-06-19 04:37:42 +0000  508)                 while 1:
a034b53b (thomie                  2010-03-19 20:09:23 +0900  509)                     newpos1 = oldpos + vector * (D1 / D12)
a034b53b (thomie                  2010-03-19 20:09:23 +0900  510)                     newpos2 = oldpos - vector * (D2 / D12)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  511)                     newpos1 = self.world.apply_boundary(newpos1)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  512)                     newpos2 = self.world.apply_boundary(newpos2)
a41e0fdc (shafi                   2008-06-19 04:37:42 +0000  513) 
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  514)                     if self.world.distance(newpos1, newpos2) >= particle_radius12:
a41e0fdc (shafi                   2008-06-19 04:37:42 +0000  515)                         break
a41e0fdc (shafi                   2008-06-19 04:37:42 +0000  516) 
a41e0fdc (shafi                   2008-06-19 04:37:42 +0000  517)                     vector *= 1.0 + 1e-7
a41e0fdc (shafi                   2008-06-19 04:37:42 +0000  518) 
f0458afe (shafi                   2007-09-26 02:45:44 +0000  519) 
f0458afe (shafi                   2007-09-26 02:45:44 +0000  520)                 # accept the new positions if there is enough space.
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  521)                 if (not self.world.check_overlap((newpos1, particle_radius1),
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  522)                                                  single.pid_particle_pair[0])) and \
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  523)                    (not self.world.check_overlap((newpos2, particle_radius2),
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  524)                                                  single.pid_particle_pair[0])):
f0458afe (shafi                   2007-09-26 02:45:44 +0000  525)                     break
f0458afe (shafi                   2007-09-26 02:45:44 +0000  526)             else:
049b6bdd (mozo                    2009-01-26 06:54:38 +0000  527)                 if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900  528)                     log.info('no space for product particles.')
4f1494cd (shafi                   2007-12-08 08:24:56 +0000  529)                 raise NoSpace()
f0458afe (shafi                   2007-09-26 02:45:44 +0000  530) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  531)             self.world.remove_particle(single.pid_particle_pair[0])
d9678fa6 (shafi                   2007-10-13 02:39:13 +0000  532) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  533)             particle1 = self.world.new_particle(product_species1.id, newpos1)
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  534)             particle2 = self.world.new_particle(product_species2.id, newpos2)
7790978e (thomie                  2010-03-19 20:27:36 +0900  535)             newsingle1 = self.create_single(particle1)
7790978e (thomie                  2010-03-19 20:27:36 +0900  536)             newsingle2 = self.create_single(particle2)
faba9096 (shafi                   2008-04-22 05:17:09 +0000  537) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  538)             self.add_single_event(newsingle1)
7790978e (thomie                  2010-03-19 20:27:36 +0900  539)             self.add_single_event(newsingle2)
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  540) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  541)             self.last_reaction = (rt, (single.pid_particle_pair[1], None),
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  542)                                   [particle1, particle2])
73377946 (shafi                   2008-08-05 23:30:41 +0000  543) 
049b6bdd (mozo                    2009-01-26 06:54:38 +0000  544)             if __debug__:
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  545)                 log.info('product1 = %s\nproduct2 = %s' % 
a034b53b (thomie                  2010-03-19 20:09:23 +0900  546)                      (str(newsingle1), str(newsingle2)))
9a84131d (shafi                   2007-09-20 02:39:40 +0000  547) 
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  548)         else:
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  549)             raise RuntimeError, 'num products >= 3 not supported.'
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  550) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  551)         self.reaction_events += 1
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  552) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  553)     def propagate_single(self, single):
3986c94d (thomie                  2010-02-05 22:29:38 +0900  554)         """The difference between a burst and a propagate is that a burst 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  555)         always takes place before the actual scheduled event for the single, 
7790978e (thomie                  2010-03-19 20:27:36 +0900  556)         while propagate_single can be called for an escape event.
3986c94d (thomie                  2010-02-05 22:29:38 +0900  557) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  558)         Another subtle difference is that burst_single always reschedules 
7790978e (thomie                  2010-03-19 20:27:36 +0900  559)         (update_event) the single, while just calling propagate does not. 
7790978e (thomie                  2010-03-19 20:27:36 +0900  560)         So whoever calls propagate_single directly should reschedule the single 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  561)         afterwards.
3986c94d (thomie                  2010-02-05 22:29:38 +0900  562) 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  563)         """
ebb84956 (moriyoshi               2009-04-15 07:44:25 +0000  564)         if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900  565)             log.debug("single.dt=%g, single.last_time=%g, self.t=%g" % (
7790978e (thomie                  2010-03-19 20:27:36 +0900  566)                single.dt, single.last_time, self.t))
c264e909 (shafi                   2008-06-23 02:14:16 +0000  567) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  568)         newpos = single.draw_new_position(single.dt, single.event_type) 
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  569)         newpos = self.world.apply_boundary(newpos)
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000  570) 
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000  571)         if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900  572)             log.debug("propagate %s: %s => %s" % (single, single.pid_particle_pair[1].position, newpos))
7393ae78 (shafi                   2008-03-20 11:47:48 +0000  573) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  574)             if self.world.check_overlap((newpos,
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  575)                                         single.pid_particle_pair[1].radius),
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  576)                                         single.pid_particle_pair[0]):
7790978e (thomie                  2010-03-19 20:27:36 +0900  577)                 raise RuntimeError('propagate_single: check_overlap failed.')
d6f015d5 (shafi                   2007-12-12 17:09:35 +0000  578) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  579)         if(single.event_type == EventType.SINGLE_REACTION and
7790978e (thomie                  2010-03-19 20:27:36 +0900  580)            single.event_type != EventType.BURST):
3986c94d (thomie                  2010-02-05 22:29:38 +0900  581)             # SINGLE_REACTION, and not a burst. No need to update, single is 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  582)             # removed anyway.
7790978e (thomie                  2010-03-19 20:27:36 +0900  583)             self.move_single_particle(single, newpos)
3986c94d (thomie                  2010-02-05 22:29:38 +0900  584)         else:
3986c94d (thomie                  2010-02-05 22:29:38 +0900  585)             # Todo. if isinstance(single, InteractionSingle):
3986c94d (thomie                  2010-02-05 22:29:38 +0900  586)             single.initialize(self.t)
7790978e (thomie                  2010-03-19 20:27:36 +0900  587)             self.move_single(single, newpos, single.pid_particle_pair[1].radius)
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  588) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  589)     def fire_single(self, single):
7790978e (thomie                  2010-03-19 20:27:36 +0900  590)         assert abs(single.dt + single.last_time - self.t) <= 1e-18 * self.t
0c22362f (shafi                   2009-12-15 05:41:14 +0000  591) 
616a0341 (thomie                  2010-02-22 17:13:13 +0900  592)         # Reaction.
7790978e (thomie                  2010-03-19 20:27:36 +0900  593)         if single.event_type == EventType.SINGLE_REACTION:
616a0341 (thomie                  2010-02-22 17:13:13 +0900  594)             if __debug__:
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  595)                 log.info('fire_single: %s' % single.event_type)
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  596)                 log.info(str(single))
0c22362f (shafi                   2009-12-15 05:41:14 +0000  597) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  598)             self.single_steps[single.event_type] += 1
8b8606aa (moriyoshi               2009-12-10 08:38:17 +0000  599) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  600)             self.propagate_single(single)
b99a9cf5 (shafi                   2007-10-31 11:45:54 +0000  601) 
f0458afe (shafi                   2007-09-26 02:45:44 +0000  602)             try:
7790978e (thomie                  2010-03-19 20:27:36 +0900  603)                 self.remove_domain(single)
7790978e (thomie                  2010-03-19 20:27:36 +0900  604)                 self.fire_single_reaction(single)
f0458afe (shafi                   2007-09-26 02:45:44 +0000  605)             except NoSpace:
b6e24d9a (thomie                  2010-02-24 15:48:36 +0900  606)                 self.reject_single_reaction(single)
b6e24d9a (thomie                  2010-02-24 15:48:36 +0900  607) 
b6e24d9a (thomie                  2010-02-24 15:48:36 +0900  608)             return
d5776c11 (shafi                   2007-09-06 23:42:04 +0000  609) 
616a0341 (thomie                  2010-02-22 17:13:13 +0900  610)         # Propagate, if not reaction.
7790978e (thomie                  2010-03-19 20:27:36 +0900  611)         single.event_type = EventType.SINGLE_ESCAPE
616a0341 (thomie                  2010-02-22 17:13:13 +0900  612)         if __debug__:
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  613)             log.info('%s: %s' % (single.event_type, str(single)))
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  614) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  615)         self.single_steps[single.event_type] += 1
616a0341 (thomie                  2010-02-22 17:13:13 +0900  616) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000  617)         # Handle immobile case first.
2f00c9c7 (shafi                   2008-07-03 03:12:19 +0000  618)         if single.getD() == 0:
a667ab59 (shafi                   2008-05-01 00:38:42 +0000  619)             # no propagation, just calculate next reaction time.
7790978e (thomie                  2010-03-19 20:27:36 +0900  620)             single.dt, single.event_type = single.determine_next_event() 
7790978e (thomie                  2010-03-19 20:27:36 +0900  621)             single.last_time = self.t
7790978e (thomie                  2010-03-19 20:27:36 +0900  622)             self.add_single_event(single)
9a603070 (shafi                   2009-11-26 06:11:41 +0000  623)             return
62e67e3d (shafi                   2007-11-30 04:45:46 +0000  624)         
3986c94d (thomie                  2010-02-05 22:29:38 +0900  625)         if single.dt != 0.0:
3986c94d (thomie                  2010-02-05 22:29:38 +0900  626)             # Propagate this particle to the exit point on the shell.
7790978e (thomie                  2010-03-19 20:27:36 +0900  627)             self.propagate_single(single)
3986c94d (thomie                  2010-02-05 22:29:38 +0900  628) 
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  629)         singlepos = single.shell.shape.position
4f1494cd (shafi                   2007-12-08 08:24:56 +0000  630) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000  631)         # (2) Clear volume.
b99a9cf5 (shafi                   2007-10-31 11:45:54 +0000  632) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  633)         min_shell = single.pid_particle_pair[1].radius * (1.0 + self.SINGLE_SHELL_FACTOR)
117beb87 (shafi                   2007-09-13 00:53:41 +0000  634) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  635)         intruders, closest, closest_distance = \
7790978e (thomie                  2010-03-19 20:27:36 +0900  636)             self.get_intruders(singlepos, min_shell, ignore=[single.domain_id, ])
3d8c644d (shafi                   2009-12-10 09:44:06 +0000  637) 
ebb84956 (moriyoshi               2009-04-15 07:44:25 +0000  638)         if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900  639)             log.debug("intruders: %s, closest: %s (dist=%g)" %\
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  640)                           (', '.join(str(i) for i in intruders), closest, closest_distance))
a667ab59 (shafi                   2008-05-01 00:38:42 +0000  641) 
3d8c644d (shafi                   2009-12-10 09:44:06 +0000  642)         if intruders:
7790978e (thomie                  2010-03-19 20:27:36 +0900  643)             burst = self.burst_non_multis(intruders)
ff8f1af9 (shafi                   2009-12-10 10:06:48 +0000  644) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  645)             obj = self.form_pair_or_multi(single, singlepos, burst)
ccdfb3a9 (shafi                   2007-09-15 02:43:33 +0000  646) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000  647)             if obj:
9a603070 (shafi                   2009-11-26 06:11:41 +0000  648)                 return
8b2a7032 (shafi                   2008-02-15 04:02:20 +0000  649) 
39da29dc (shafi                   2008-04-25 06:02:00 +0000  650)             # if nothing was formed, recheck closest and restore shells.
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  651)             burst = uniq(burst)
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  652) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  653)             closest, closest_distance = \
7790978e (thomie                  2010-03-19 20:27:36 +0900  654)                 self.get_closest_obj(singlepos, ignore = [single.domain_id, ])
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  655)             self.update_single(single, closest, closest_distance)
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  656)             for s in burst:
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  657)                 if not isinstance(s, Single):
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  658)                     continue
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  659)                 assert s.is_reset()
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  660)                 closest, closest_distance = self.get_closest_obj(
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  661)                     s.shell.shape.position, ignore = [s.domain_id, ])
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  662) 
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  663)                 self.update_single(s, closest, closest_distance)
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  664)                 self.update_single_event(self.t + s.dt, s)
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  665)                 if __debug__:
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  666)                     log.debug('update single: %s\nclosest %s, closest_distance %g' %
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  667)                           (repr(s), closest, closest_distance))
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  668)         else:
cd391fac (Moriyoshi Koizumi       2010-06-25 16:30:50 +0900  669)             self.update_single(single, closest, closest_distance)
117beb87 (shafi                   2007-09-13 00:53:41 +0000  670) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  671)         self.add_single_event(single)
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  672) 
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  673)         if __debug__:
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  674)             log.debug('update shell: %s.' % repr(single))
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  675)             
9a603070 (shafi                   2009-11-26 06:11:41 +0000  676)         return
b99a9cf5 (shafi                   2007-10-31 11:45:54 +0000  677) 
b6e24d9a (thomie                  2010-02-24 15:48:36 +0900  678)     def reject_single_reaction(self, single):
b6e24d9a (thomie                  2010-02-24 15:48:36 +0900  679)         if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900  680)             log.info('single reaction; placing product failed.')
b6e24d9a (thomie                  2010-02-24 15:48:36 +0900  681)         self.domains[single.domain_id] = single
7790978e (thomie                  2010-03-19 20:27:36 +0900  682)         self.move_shell(single.shell_id_shell_pair)
7790978e (thomie                  2010-03-19 20:27:36 +0900  683)         self.rejected_moves += 1
987c4137 (thomie                  2010-03-01 13:32:46 +0900  684)         single.initialize(self.t)
7790978e (thomie                  2010-03-19 20:27:36 +0900  685)         self.add_single_event(single)
b6e24d9a (thomie                  2010-02-24 15:48:36 +0900  686) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  687)     def calculate_single_shell_size(self, single, closest, 
7790978e (thomie                  2010-03-19 20:27:36 +0900  688)                                  distance, shell_distance):
a034b53b (thomie                  2010-03-19 20:09:23 +0900  689)         assert isinstance(closest, Single)
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  690) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  691)         min_radius1 = single.pid_particle_pair[1].radius
68760ef5 (shafi                   2008-06-22 02:16:31 +0000  692)         D1 = single.getD()
68760ef5 (shafi                   2008-06-22 02:16:31 +0000  693) 
68760ef5 (shafi                   2008-06-22 02:16:31 +0000  694)         if D1 == 0:
7790978e (thomie                  2010-03-19 20:27:36 +0900  695)             return min_radius1
68760ef5 (shafi                   2008-06-22 02:16:31 +0000  696) 
68760ef5 (shafi                   2008-06-22 02:16:31 +0000  697)         D2 = closest.getD()
7790978e (thomie                  2010-03-19 20:27:36 +0900  698)         min_radius2 = closest.pid_particle_pair[1].radius
7790978e (thomie                  2010-03-19 20:27:36 +0900  699)         min_radius12 = min_radius1 + min_radius2
a034b53b (thomie                  2010-03-19 20:09:23 +0900  700)         sqrtD1 = math.sqrt(D1)
68760ef5 (shafi                   2008-06-22 02:16:31 +0000  701)             
7790978e (thomie                  2010-03-19 20:27:36 +0900  702)         shell_size = min(sqrtD1 / (sqrtD1 + math.sqrt(D2))
7790978e (thomie                  2010-03-19 20:27:36 +0900  703)                         * (distance - min_radius12) + min_radius1,
7790978e (thomie                  2010-03-19 20:27:36 +0900  704)                         shell_distance / SAFETY)
7790978e (thomie                  2010-03-19 20:27:36 +0900  705)         if shell_size < min_radius1:
7790978e (thomie                  2010-03-19 20:27:36 +0900  706)             shell_size = min_radius1
68760ef5 (shafi                   2008-06-22 02:16:31 +0000  707) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  708)         return shell_size
68760ef5 (shafi                   2008-06-22 02:16:31 +0000  709) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  710)     def update_single(self, single, closest, distance_to_shell): 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  711)         # Todo. assert not isinstance(single, InteractionSingle)
3986c94d (thomie                  2010-02-05 22:29:38 +0900  712) 
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  713)         singlepos = single.shell.shape.position
a034b53b (thomie                  2010-03-19 20:09:23 +0900  714)         if isinstance(closest, Single):
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  715)             closestpos = closest.shell.shape.position
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  716)             distance_to_closest = self.world.distance(singlepos, closestpos)
7790978e (thomie                  2010-03-19 20:27:36 +0900  717)             new_shell_size = self.calculate_single_shell_size(single, closest, 
7790978e (thomie                  2010-03-19 20:27:36 +0900  718)                                                       distance_to_closest,
7790978e (thomie                  2010-03-19 20:27:36 +0900  719)                                                       distance_to_shell)
faba9096 (shafi                   2008-04-22 05:17:09 +0000  720)         else:  # Pair or Multi
7790978e (thomie                  2010-03-19 20:27:36 +0900  721)             new_shell_size = distance_to_shell / SAFETY
a034b53b (thomie                  2010-03-19 20:09:23 +0900  722)             new_shell_size = max(new_shell_size, single.pid_particle_pair[1].radius)
a7e32062 (shafi                   2007-09-08 02:27:45 +0000  723) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  724)         new_shell_size = min(new_shell_size, self.get_max_shell_size())
e3e19d34 (shafi                   2007-11-28 23:28:03 +0000  725) 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  726)         # Resize shell, don't change position.
7790978e (thomie                  2010-03-19 20:27:36 +0900  727)         # Note: this should be done before determine_next_event.
7790978e (thomie                  2010-03-19 20:27:36 +0900  728)         self.move_single_shell(single, singlepos, new_shell_size)        
3986c94d (thomie                  2010-02-05 22:29:38 +0900  729) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  730)         single.dt, single.event_type = single.determine_next_event()
7790978e (thomie                  2010-03-19 20:27:36 +0900  731)         single.last_time = self.t
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  732) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  733)     def fire_pair(self, pair):
7790978e (thomie                  2010-03-19 20:27:36 +0900  734)         assert self.check_obj(pair)
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  735) 
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  736)         single1 = pair.single1
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  737)         single2 = pair.single2
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  738)         particle1 = single1.pid_particle_pair
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  739)         particle2 = single2.pid_particle_pair
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  740)         pos1 = particle1[1].position
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  741)         pos2 = particle2[1].position
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  742)         
7790978e (thomie                  2010-03-19 20:27:36 +0900  743)         if pair.event_type == EventType.IV_EVENT:
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  744)             # Draw actual pair event for iv at very last minute.
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  745)             r0 = self.world.distance(pos1, pos2)
7790978e (thomie                  2010-03-19 20:27:36 +0900  746)             pair.event_type = pair.draw_iv_event_type(r0)
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  747) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  748)         self.pair_steps[pair.event_type] += 1
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  749) 
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  750)         if __debug__:
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  751)             log.info('%s: %s, %s' % (pair.event_type, single1, single2))
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  752) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  753)         old_com = pair.com
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  754) 
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  755)         # Four cases:
7e3d1fd5 (thomie                  2010-02-23 15:27:07 +0900  756)         #  1. Single reaction
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  757)         #  2. Pair reaction
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  758)         #  3a. IV escape
7790978e (thomie                  2010-03-19 20:27:36 +0900  759)         #  3b. com escape
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  760) 
7e3d1fd5 (thomie                  2010-02-23 15:27:07 +0900  761)         #
7e3d1fd5 (thomie                  2010-02-23 15:27:07 +0900  762)         # 1. Single reaction
7e3d1fd5 (thomie                  2010-02-23 15:27:07 +0900  763)         #
7790978e (thomie                  2010-03-19 20:27:36 +0900  764)         if pair.event_type == EventType.SINGLE_REACTION:
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  765)             reactingsingle = pair.reactingsingle
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  766) 
049b6bdd (mozo                    2009-01-26 06:54:38 +0000  767)             if __debug__:
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  768)                 log.info('reactant = %s' % str(reactingsingle))
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  769) 
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  770)             if reactingsingle == single1:
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  771)                 theothersingle = single2
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  772)             else:
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  773)                 theothersingle = single1
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  774) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  775)             self.burst_pair(pair)
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  776) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  777)             self.add_single_event(theothersingle)
3c3e55d8 (shafi                   2008-02-19 00:44:24 +0000  778) 
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  779)             try:
7790978e (thomie                  2010-03-19 20:27:36 +0900  780)                 self.remove_domain(reactingsingle)
7790978e (thomie                  2010-03-19 20:27:36 +0900  781)                 self.fire_single_reaction(reactingsingle)
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  782)             except NoSpace:
b6e24d9a (thomie                  2010-02-24 15:48:36 +0900  783)                 self.reject_single_reaction(reactingsingle)
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  784) 
9a603070 (shafi                   2009-11-26 06:11:41 +0000  785)             return
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  786)         
7e3d1fd5 (thomie                  2010-02-23 15:27:07 +0900  787)         #
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  788)         # 2. Pair reaction
7e3d1fd5 (thomie                  2010-02-23 15:27:07 +0900  789)         #
7790978e (thomie                  2010-03-19 20:27:36 +0900  790)         if pair.event_type == EventType.IV_REACTION:
a034b53b (thomie                  2010-03-19 20:09:23 +0900  791)             if len(pair.rt.products) == 1:
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  792)                 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  793)                 species3 = self.world.get_species(pair.rt.products[0])
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  794) 
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  795)                 # calculate new R
7790978e (thomie                  2010-03-19 20:27:36 +0900  796)                 event_type = pair.event_type
7790978e (thomie                  2010-03-19 20:27:36 +0900  797)                 new_com = pair.draw_new_com(pair.dt, event_type)
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  798)                 
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  799)                 if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900  800)                     shell_size = pair.get_shell_size()
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  801)                     assert self.world.distance(old_com, new_com) + species3.radius <\
7790978e (thomie                  2010-03-19 20:27:36 +0900  802)                         shell_size
d6f015d5 (shafi                   2007-12-12 17:09:35 +0000  803) 
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  804)                 new_com = self.world.apply_boundary(new_com)
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  805) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  806)                 self.world.remove_particle(single1.pid_particle_pair[0])
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  807)                 self.world.remove_particle(single2.pid_particle_pair[0])
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  808) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  809)                 particle = self.world.new_particle(species3.id, new_com)
7790978e (thomie                  2010-03-19 20:27:36 +0900  810)                 newsingle = self.create_single(particle)
7790978e (thomie                  2010-03-19 20:27:36 +0900  811)                 self.add_single_event(newsingle)
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  812) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  813)                 self.reaction_events += 1
73377946 (shafi                   2008-08-05 23:30:41 +0000  814) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  815)                 self.last_reaction = (pair.rt, (particle1, particle2),
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  816)                                       [particle])
73377946 (shafi                   2008-08-05 23:30:41 +0000  817) 
049b6bdd (mozo                    2009-01-26 06:54:38 +0000  818)                 if __debug__:
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  819)                     log.info('product = %s' % str(newsingle))
73377946 (shafi                   2008-08-05 23:30:41 +0000  820) 
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  821)             else:
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  822)                 raise NotImplementedError,\
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  823)                       'num products >= 2 not supported.'
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  824) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  825)             self.remove_domain(pair)
d6f015d5 (shafi                   2007-12-12 17:09:35 +0000  826) 
9a603070 (shafi                   2009-11-26 06:11:41 +0000  827)             return
81c0bdc7 (shafi                   2007-09-22 02:26:53 +0000  828) 
81c0bdc7 (shafi                   2007-09-22 02:26:53 +0000  829)         #
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  830)         # 3a. Escaping through a_r.
91a21f9f (thomie                  2010-02-25 20:49:52 +0900  831)         # 3b. Escaping through a_R.
81c0bdc7 (shafi                   2007-09-22 02:26:53 +0000  832)         #
7790978e (thomie                  2010-03-19 20:27:36 +0900  833)         elif(pair.event_type == EventType.IV_ESCAPE or
7790978e (thomie                  2010-03-19 20:27:36 +0900  834)              pair.event_type == EventType.COM_ESCAPE):
3a467150 (thomie                  2010-02-23 17:43:07 +0900  835)             dt = pair.dt
7790978e (thomie                  2010-03-19 20:27:36 +0900  836)             event_type = pair.event_type
7790978e (thomie                  2010-03-19 20:27:36 +0900  837)             single1, single2 = self.propagate_pair(pair, dt, event_type)
7790978e (thomie                  2010-03-19 20:27:36 +0900  838)             self.add_single_event(single1)
7790978e (thomie                  2010-03-19 20:27:36 +0900  839)             self.add_single_event(single2)
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  840)         else:
7790978e (thomie                  2010-03-19 20:27:36 +0900  841)             raise SystemError, 'Bug: invalid event_type.'
c06b46b5 (shafi                   2007-06-29 02:21:13 +0000  842) 
9a603070 (shafi                   2009-11-26 06:11:41 +0000  843)         return
b99a9cf5 (shafi                   2007-10-31 11:45:54 +0000  844) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  845)     def fire_multi(self, multi):
f86d6843 (shafi                   2009-12-21 07:21:26 +0000  846)         self.multi_steps[2] += 1  # multi_steps[2]: total multi steps
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  847)         multi.step()
0c22362f (shafi                   2009-12-15 05:41:14 +0000  848) 
8b8606aa (moriyoshi               2009-12-10 08:38:17 +0000  849)         if __debug__:
00000000 (Not Committed Yet       2010-08-10 13:37:17 +0200  850)             log.info(multi.last_event)
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  851) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  852)         if multi.last_event == EventType.MULTI_REACTION:
7790978e (thomie                  2010-03-19 20:27:36 +0900  853)             self.reaction_events += 1
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  854)             self.last_reaction = multi.last_reaction
cb68c704 (shafi                   2008-07-18 20:26:51 +0000  855) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  856)         if multi.last_event is not None:
7790978e (thomie                  2010-03-19 20:27:36 +0900  857)             self.break_up_multi(multi)
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  858)             self.multi_steps[multi.last_event] += 1
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  859)         else:
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  860)             self.add_multi_event(multi)
9a603070 (shafi                   2009-11-26 06:11:41 +0000  861) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  862)     def break_up_multi(self, multi):
7790978e (thomie                  2010-03-19 20:27:36 +0900  863)         self.remove_domain(multi)
faba9096 (shafi                   2008-04-22 05:17:09 +0000  864) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000  865)         singles = []
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  866)         for pid_particle_pair in multi.particles:
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  867)             single = self.create_single(pid_particle_pair)
7790978e (thomie                  2010-03-19 20:27:36 +0900  868)             self.add_single_event(single)
a034b53b (thomie                  2010-03-19 20:09:23 +0900  869)             singles.append(single)
faba9096 (shafi                   2008-04-22 05:17:09 +0000  870) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000  871)         return singles
3702a64e (shafi                   2008-03-28 00:26:28 +0000  872) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  873)     def burst_multi(self, multi):
b9cb9048 (shafi                   2008-07-22 23:22:35 +0000  874)         #multi.sim.sync()
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000  875)         assert isinstance(multi, Multi)
7790978e (thomie                  2010-03-19 20:27:36 +0900  876)         singles = self.break_up_multi(multi)
b76b4839 (shafi                   2008-03-27 06:52:00 +0000  877) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000  878)         return singles
b76b4839 (shafi                   2008-03-27 06:52:00 +0000  879) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  880)     def burst_single(self, single):
7790978e (thomie                  2010-03-19 20:27:36 +0900  881)         assert self.t >= single.last_time
7790978e (thomie                  2010-03-19 20:27:36 +0900  882)         assert self.t <= single.last_time + single.dt
4f1494cd (shafi                   2007-12-08 08:24:56 +0000  883) 
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  884)         oldpos = single.shell.shape.position
bf8ee670 (thomie                  2010-02-26 22:08:33 +0900  885)         old_shell_size = single.get_shell_size()
bf8ee670 (thomie                  2010-02-26 22:08:33 +0900  886) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  887)         particle_radius = single.pid_particle_pair[1].radius
7393ae78 (shafi                   2008-03-20 11:47:48 +0000  888) 
8d6e2c81 (thomie                  2010-02-24 00:40:18 +0900  889)         # Override dt, burst happens before single's scheduled event.
7790978e (thomie                  2010-03-19 20:27:36 +0900  890)         single.dt = self.t - single.last_time
7790978e (thomie                  2010-03-19 20:27:36 +0900  891)         # Override event_type. Always call gf.drawR on BURST.
7790978e (thomie                  2010-03-19 20:27:36 +0900  892)         single.event_type = EventType.BURST
7790978e (thomie                  2010-03-19 20:27:36 +0900  893)         self.propagate_single(single)
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  894) 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  895)         newpos = single.pid_particle_pair[1].position
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  896)         assert self.world.distance(newpos, oldpos) <= old_shell_size - particle_radius
7790978e (thomie                  2010-03-19 20:27:36 +0900  897)         # Displacement check is in NonInteractionSingle.draw_new_position.
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  898) 
3986c94d (thomie                  2010-02-05 22:29:38 +0900  899)         # Todo. if isinstance(single, InteractionSingle):
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900  900)         self.update_single_event(self.t, single)
18f25746 (shafi                   2007-03-29 03:44:36 +0000  901) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  902)         assert single.shell.shape.radius == particle_radius
5be27630 (thomie                  2010-02-05 14:36:48 +0900  903) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  904)     def burst_pair(self, pair):
e4d9fa0c (shafi@senin             2010-01-14 22:40:25 +0900  905)         if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900  906)             log.debug('burst_pair: %s', pair)
e4d9fa0c (shafi@senin             2010-01-14 22:40:25 +0900  907) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  908)         assert self.t >= pair.last_time
7790978e (thomie                  2010-03-19 20:27:36 +0900  909)         assert self.t <= pair.last_time + pair.dt
4f1494cd (shafi                   2007-12-08 08:24:56 +0000  910) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  911)         dt = self.t - pair.last_time 
7790978e (thomie                  2010-03-19 20:27:36 +0900  912)         # Override event_type. Always call sgf.drawR and pgf.drawR on BURST.
7790978e (thomie                  2010-03-19 20:27:36 +0900  913)         event_type = EventType.BURST
7790978e (thomie                  2010-03-19 20:27:36 +0900  914)         single1, single2 = self.propagate_pair(pair, dt, event_type)
4c0efc88 (thomie                  2010-02-23 22:39:15 +0900  915) 
4c0efc88 (thomie                  2010-02-23 22:39:15 +0900  916)         return single1, single2
4c0efc88 (thomie                  2010-02-23 22:39:15 +0900  917) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  918)     def propagate_pair(self, pair, dt, event_type):
e4d9fa0c (shafi@senin             2010-01-14 22:40:25 +0900  919)         single1 = pair.single1
e4d9fa0c (shafi@senin             2010-01-14 22:40:25 +0900  920)         single2 = pair.single2
e4d9fa0c (shafi@senin             2010-01-14 22:40:25 +0900  921) 
5be27630 (thomie                  2010-02-05 14:36:48 +0900  922)         particle1 = single1.pid_particle_pair
5be27630 (thomie                  2010-02-05 14:36:48 +0900  923)         particle2 = single2.pid_particle_pair
4f1494cd (shafi                   2007-12-08 08:24:56 +0000  924) 
3cd4a83c (thomie                  2010-02-23 23:01:36 +0900  925)         pos1 = particle1[1].position
3cd4a83c (thomie                  2010-02-23 23:01:36 +0900  926)         pos2 = particle2[1].position
4f1494cd (shafi                   2007-12-08 08:24:56 +0000  927) 
5be27630 (thomie                  2010-02-05 14:36:48 +0900  928)         if dt > 0.0:
210d665f (shafi                   2009-12-15 09:08:50 +0000  929)             D1 = particle1[1].D
b69f926a (shafi                   2009-12-15 09:30:43 +0000  930)             D2 = particle2[1].D
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000  931) 
acd00e2c (Moriyoshi Koizumi       2010-03-11 14:21:35 +0900  932)             pos2t = self.world.cyclic_transpose(pos2, pos1)
7790978e (thomie                  2010-03-19 20:27:36 +0900  933)             old_inter_particle = pos2t - pos1
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  934)             r0 = self.world.distance(pos1, pos2)
7790978e (thomie                  2010-03-19 20:27:36 +0900  935)             assert feq(r0, length(old_inter_particle))
1c710015 (shafi                   2009-06-26 04:47:18 +0000  936) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  937)             old_com = pair.com
1c710015 (shafi                   2009-06-26 04:47:18 +0000  938) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  939)             newpos1, newpos2 = pair.draw_new_positions(dt, r0, 
7790978e (thomie                  2010-03-19 20:27:36 +0900  940)                                                      old_inter_particle, 
7790978e (thomie                  2010-03-19 20:27:36 +0900  941)                                                      event_type)
17e1fe68 (shafi@tcs1.gsc.riken.jp 2010-02-03 12:26:26 +0900  942) 
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  943)             newpos1 = self.world.apply_boundary(newpos1)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900  944)             newpos2 = self.world.apply_boundary(newpos2)
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  945)             assert not self.world.check_overlap((newpos1, particle1[1].radius),
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  946)                                                 particle1[0], particle2[0])
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  947)             assert not self.world.check_overlap((newpos2, particle2[1].radius),
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900  948)                                                 particle1[0], particle2[0])
7790978e (thomie                  2010-03-19 20:27:36 +0900  949)             assert self.check_pair_pos(pair, newpos1, newpos2, old_com,\
bf8ee670 (thomie                  2010-02-26 22:08:33 +0900  950)                                          pair.get_shell_size())
37a895c1 (Moriyoshi Koizumi       2010-02-08 19:40:00 +0900  951)         else:
37a895c1 (Moriyoshi Koizumi       2010-02-08 19:40:00 +0900  952)             newpos1 = particle1[1].position
37a895c1 (Moriyoshi Koizumi       2010-02-08 19:40:00 +0900  953)             newpos2 = particle2[1].position
4f1494cd (shafi                   2007-12-08 08:24:56 +0000  954) 
3cd4a83c (thomie                  2010-02-23 23:01:36 +0900  955)         if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900  956)             log.debug("fire_pair: #1 { %s: %s => %s }" % (single1, pos1, newpos1))
7790978e (thomie                  2010-03-19 20:27:36 +0900  957)             log.debug("fire_pair: #2 { %s: %s => %s }" % (single2, pos2, newpos2))
3cd4a83c (thomie                  2010-02-23 23:01:36 +0900  958) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900  959)         single1.initialize(self.t)
a034b53b (thomie                  2010-03-19 20:09:23 +0900  960)         single2.initialize(self.t)
df083205 (shafi                   2007-09-19 23:30:15 +0000  961)         
7790978e (thomie                  2010-03-19 20:27:36 +0900  962)         self.remove_domain(pair)
7edbf48b (moriyoshi               2009-12-22 06:34:47 +0000  963)         assert single1.domain_id not in self.domains
7edbf48b (moriyoshi               2009-12-22 06:34:47 +0000  964)         assert single2.domain_id not in self.domains
7edbf48b (moriyoshi               2009-12-22 06:34:47 +0000  965)         self.domains[single1.domain_id] = single1
7edbf48b (moriyoshi               2009-12-22 06:34:47 +0000  966)         self.domains[single2.domain_id] = single2
7790978e (thomie                  2010-03-19 20:27:36 +0900  967)         self.move_single(single1, newpos1, particle1[1].radius)
7790978e (thomie                  2010-03-19 20:27:36 +0900  968)         self.move_single(single2, newpos2, particle2[1].radius)
d6f015d5 (shafi                   2007-12-12 17:09:35 +0000  969) 
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900  970)         if __debug__:
e9d894b2 (thomie                  2010-02-26 21:42:07 +0900  971)             container = self.get_container(single1.shell)
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  972)             assert container[single1.shell_id].shape.radius == single1.shell.shape.radius
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  973)             assert container[single2.shell_id].shape.radius == single2.shell.shape.radius
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900  974) 
737c6f80 (thomie                  2010-03-13 00:19:25 +0900  975)             if type(single1.shell) is CylindricalShell:
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  976)                 assert container[single1.shell_id].shape.size == single1.shell.shape.size
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  977)                 assert container[single2.shell_id].shape.size == single2.shell.shape.size
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900  978) 
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  979)         assert single1.shell.shape.radius == particle1[1].radius
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900  980)         assert single2.shell.shape.radius == particle2[1].radius
e4d9fa0c (shafi@senin             2010-01-14 22:40:25 +0900  981) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  982)         assert self.check_obj(single1)
7790978e (thomie                  2010-03-19 20:27:36 +0900  983)         assert self.check_obj(single2)
3cd4a83c (thomie                  2010-02-23 23:01:36 +0900  984) 
d0cf1c00 (shafi                   2007-11-30 02:43:37 +0000  985)         return single1, single2
8be9229e (shafi                   2007-03-21 03:35:34 +0000  986) 
7790978e (thomie                  2010-03-19 20:27:36 +0900  987)     def form_pair_or_multi(self, single, singlepos, neighbors):
faba9096 (shafi                   2008-04-22 05:17:09 +0000  988)         assert neighbors
fd242469 (shafi                   2007-07-04 02:42:21 +0000  989) 
ce668d73 (shafi                   2009-12-18 07:42:24 +0000  990)         # sort burst neighbors by distance
7790978e (thomie                  2010-03-19 20:27:36 +0900  991)         dists = self.obj_distance_array(singlepos, neighbors)
ead98343 (shafi                   2009-12-18 09:32:03 +0000  992)         if len(dists) >= 2:
ead98343 (shafi                   2009-12-18 09:32:03 +0000  993)             n = dists.argsort()
ead98343 (shafi                   2009-12-18 09:32:03 +0000  994)             dists = dists.take(n)
ead98343 (shafi                   2009-12-18 09:32:03 +0000  995)             neighbors = numpy.take(neighbors, n)
ce668d73 (shafi                   2009-12-18 07:42:24 +0000  996) 
ce668d73 (shafi                   2009-12-18 07:42:24 +0000  997)         # First, try forming a Pair.
ce668d73 (shafi                   2009-12-18 07:42:24 +0000  998)         if isinstance(neighbors[0], Single):
7790978e (thomie                  2010-03-19 20:27:36 +0900  999)             obj = self.form_pair(single, singlepos, neighbors[0], neighbors[1:])
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1000)             if obj:
bc0d5aa9 (shafi                   2009-05-07 07:46:16 +0000 1001)                 return obj
117beb87 (shafi                   2007-09-13 00:53:41 +0000 1002) 
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1003)         # If a Pair is not formed, then try forming a Multi.
7790978e (thomie                  2010-03-19 20:27:36 +0900 1004)         obj = self.form_multi(single, neighbors, dists)
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1005)         if obj:
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1006)             return obj
d3d23575 (shafi                   2008-04-23 02:39:49 +0000 1007) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1008) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1009)     def form_pair(self, single1, pos1, single2, burst):
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1010)         if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1011)            log.debug('trying to form %s' %
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1012)                  'Pair(%s, %s)' % (single1.pid_particle_pair, 
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1013)                                      single2.pid_particle_pair))
776de2f2 (shafi                   2008-05-04 23:29:34 +0000 1014) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1015)         assert single1.is_reset()
7790978e (thomie                  2010-03-19 20:27:36 +0900 1016)         assert single2.is_reset()
ec1b4735 (shafi                   2007-04-26 04:20:06 +0000 1017) 
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000 1018)         radius1 = single1.pid_particle_pair[1].radius
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000 1019)         radius2 = single2.pid_particle_pair[1].radius
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1020) 
9afaf84f (shafi                   2008-05-06 04:30:58 +0000 1021)         sigma = radius1 + radius2
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1022) 
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000 1023)         D1, D2 = single1.pid_particle_pair[1].D, single2.pid_particle_pair[1].D
75e8dfa0 (shafi                   2007-06-30 04:22:03 +0000 1024)         D12 = D1 + D2
36dbd52e (shafi                   2007-08-27 16:49:29 +0000 1025) 
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900 1026)         assert (pos1 - single1.shell.shape.position).sum() == 0
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900 1027)         pos2 = single2.shell.shape.position
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900 1028)         r0 = self.world.distance(pos1, pos2)
7790978e (thomie                  2010-03-19 20:27:36 +0900 1029)         distance_from_sigma = r0 - sigma
7790978e (thomie                  2010-03-19 20:27:36 +0900 1030)         assert distance_from_sigma >= 0, \
7790978e (thomie                  2010-03-19 20:27:36 +0900 1031)             ('distance_from_sigma (pair gap) between %s and %s = %g < 0' %
7790978e (thomie                  2010-03-19 20:27:36 +0900 1032)              (single1, single2, distance_from_sigma))
7790978e (thomie                  2010-03-19 20:27:36 +0900 1033) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1034)         shell_size1 = r0 * D1 / D12 + radius1
7790978e (thomie                  2010-03-19 20:27:36 +0900 1035)         shell_size2 = r0 * D2 / D12 + radius2
7790978e (thomie                  2010-03-19 20:27:36 +0900 1036)         shell_size_margin1 = radius1 * 2 #* self.SINGLE_SHELL_FACTOR
7790978e (thomie                  2010-03-19 20:27:36 +0900 1037)         shell_size_margin2 = radius2 * 2 #* self.SINGLE_SHELL_FACTOR
7790978e (thomie                  2010-03-19 20:27:36 +0900 1038)         shell_size_with_margin1 = shell_size1 + shell_size_margin1
7790978e (thomie                  2010-03-19 20:27:36 +0900 1039)         shell_size_with_margin2 = shell_size2 + shell_size_margin2
7790978e (thomie                  2010-03-19 20:27:36 +0900 1040)         if shell_size_with_margin1  >= shell_size_with_margin2:
7790978e (thomie                  2010-03-19 20:27:36 +0900 1041)             min_shell_size = shell_size1
7790978e (thomie                  2010-03-19 20:27:36 +0900 1042)             shell_size_margin = shell_size_margin1
9afaf84f (shafi                   2008-05-06 04:30:58 +0000 1043)         else:
7790978e (thomie                  2010-03-19 20:27:36 +0900 1044)             min_shell_size = shell_size2
7790978e (thomie                  2010-03-19 20:27:36 +0900 1045)             shell_size_margin = shell_size_margin2
9487b2a7 (shafi                   2008-04-25 01:50:16 +0000 1046) 
9487b2a7 (shafi                   2008-04-25 01:50:16 +0000 1047)         # 1. Shell cannot be larger than max shell size or sim cell size.
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1048)         com = self.world.calculate_pair_CoM(pos1, pos2, D1, D2)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900 1049)         com = self.world.apply_boundary(com)
7790978e (thomie                  2010-03-19 20:27:36 +0900 1050)         min_shell_size_with_margin = min_shell_size + shell_size_margin
7790978e (thomie                  2010-03-19 20:27:36 +0900 1051)         max_shell_size = min(self.get_max_shell_size(),
7790978e (thomie                  2010-03-19 20:27:36 +0900 1052)                            distance_from_sigma * 100 + sigma + shell_size_margin)
84b7b61b (shafi                   2007-04-26 02:05:20 +0000 1053) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1054)         if min_shell_size_with_margin >= max_shell_size:
049b6bdd (mozo                    2009-01-26 06:54:38 +0000 1055)             if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900 1056)                 log.debug('%s not formed: min_shell_size %g >= max_shell_size %g' %
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1057)                           ('Pair(%s, %s)' % (single1.pid_particle_pair[0], 
c9f960e9 (thomie                  2010-02-27 03:01:51 +0900 1058)                                                single2.pid_particle_pair[0]),
7790978e (thomie                  2010-03-19 20:27:36 +0900 1059)                            min_shell_size_with_margin, max_shell_size))
297f3a65 (shafi                   2008-05-01 02:19:07 +0000 1060)             return None
297f3a65 (shafi                   2008-05-01 02:19:07 +0000 1061) 
1c710015 (shafi                   2009-06-26 04:47:18 +0000 1062)         # Here, we have to take into account of the burst Singles in
297f3a65 (shafi                   2008-05-01 02:19:07 +0000 1063)         # this step.  The simple check for closest below could miss
297f3a65 (shafi                   2008-05-01 02:19:07 +0000 1064)         # some of them, because sizes of these Singles for this
297f3a65 (shafi                   2008-05-01 02:19:07 +0000 1065)         # distance check has to include SINGLE_SHELL_FACTOR, while
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1066)         # these burst objects have zero mobility radii.  This is not
297f3a65 (shafi                   2008-05-01 02:19:07 +0000 1067)         # beautiful, a cleaner framework may be possible.
f356bdeb (shafi                   2008-04-30 01:53:48 +0000 1068) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1069)         closest, closest_shell_distance = None, numpy.inf
ff8f1af9 (shafi                   2009-12-10 10:06:48 +0000 1070)         for b in burst:
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1071)             if isinstance(b, Single):
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900 1072)                 bpos = b.shell.shape.position
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900 1073)                 d = self.world.distance(com, bpos) \
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1074)                     - b.pid_particle_pair[1].radius * (1.0 + self.SINGLE_SHELL_FACTOR)
7790978e (thomie                  2010-03-19 20:27:36 +0900 1075)                 if d < closest_shell_distance:
7790978e (thomie                  2010-03-19 20:27:36 +0900 1076)                     closest, closest_shell_distance = b, d
1cff6ee8 (shafi                   2008-05-01 02:34:29 +0000 1077) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1078)         if closest_shell_distance <= min_shell_size_with_margin:
049b6bdd (mozo                    2009-01-26 06:54:38 +0000 1079)             if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1080)                 log.debug('%s not formed: squeezed by burst neighbor %s' %
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1081)                       ('Pair(%s, %s)' % (single1.pid_particle_pair[0], 
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1082)                                            single2.pid_particle_pair[0]),
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1083)                        closest))
1cff6ee8 (shafi                   2008-05-01 02:34:29 +0000 1084)             return None
a667ab59 (shafi                   2008-05-01 00:38:42 +0000 1085) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1086)         assert closest_shell_distance > 0
7790978e (thomie                  2010-03-19 20:27:36 +0900 1087)         c, d = self.get_closest_obj(com, ignore=[single1.domain_id, single2.domain_id])
7790978e (thomie                  2010-03-19 20:27:36 +0900 1088)         if d < closest_shell_distance:
7790978e (thomie                  2010-03-19 20:27:36 +0900 1089)             closest, closest_shell_distance = c, d
f356bdeb (shafi                   2008-04-30 01:53:48 +0000 1090) 
049b6bdd (mozo                    2009-01-26 06:54:38 +0000 1091)         if __debug__:
7790978e (thomie                  2010-03-19 20:27:36 +0900 1092)             log.debug('Pair closest neighbor: %s %g, min_shell_with_margin %g' %
7790978e (thomie                  2010-03-19 20:27:36 +0900 1093)                   (closest, closest_shell_distance, min_shell_size_with_margin))
f356bdeb (shafi                   2008-04-30 01:53:48 +0000 1094) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1095)         assert closest_shell_distance > 0
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1096) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1097)         if isinstance(closest, Single):
297f3a65 (shafi                   2008-05-01 02:19:07 +0000 1098) 
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000 1099)             D_closest = closest.pid_particle_pair[1].D
297f3a65 (shafi                   2008-05-01 02:19:07 +0000 1100)             D_tot = D_closest + D12
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900 1101)             closest_distance = self.world.distance(com, closest.pid_particle_pair[1].position) ##??
297f3a65 (shafi                   2008-05-01 02:19:07 +0000 1102) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1103)             closest_min_radius = closest.pid_particle_pair[1].radius
7790978e (thomie                  2010-03-19 20:27:36 +0900 1104)             closest_min_shell = closest_min_radius * \
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1105)                 (self.SINGLE_SHELL_FACTOR + 1.0)
c814556f (shafi                   2008-04-27 01:47:29 +0000 1106) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1107)             shell_size = min((D12 / D_tot) *
7790978e (thomie                  2010-03-19 20:27:36 +0900 1108)                             (closest_distance - min_shell_size 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1109)                              - closest_min_radius) + min_shell_size,
7790978e (thomie                  2010-03-19 20:27:36 +0900 1110)                             closest_distance - closest_min_shell,
7790978e (thomie                  2010-03-19 20:27:36 +0900 1111)                             closest_shell_distance)
297f3a65 (shafi                   2008-05-01 02:19:07 +0000 1112) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1113)             shell_size /= SAFETY
7790978e (thomie                  2010-03-19 20:27:36 +0900 1114)             assert shell_size < closest_shell_distance
c814556f (shafi                   2008-04-27 01:47:29 +0000 1115) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1116)         else:
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1117)             assert isinstance(closest, (Pair, Multi, None.__class__))
c814556f (shafi                   2008-04-27 01:47:29 +0000 1118) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1119)             shell_size = closest_shell_distance / SAFETY
297f3a65 (shafi                   2008-05-01 02:19:07 +0000 1120) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1121)         if shell_size <= min_shell_size_with_margin:
049b6bdd (mozo                    2009-01-26 06:54:38 +0000 1122)             if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1123)                 log.debug('%s not formed: squeezed by %s' %
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1124)                       ('Pair(%s, %s)' % (single1.pid_particle_pair[0], 
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1125)                                            single2.pid_particle_pair[0]),
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1126)                        closest))
7a65bbc4 (shafi                   2008-05-01 03:12:58 +0000 1127)             return None
7a65bbc4 (shafi                   2008-05-01 03:12:58 +0000 1128) 
7a65bbc4 (shafi                   2008-05-01 03:12:58 +0000 1129) 
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900 1130)         d1 = self.world.distance(com, pos1)
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900 1131)         d2 = self.world.distance(com, pos2)
7b98afc5 (shafi                   2008-04-29 06:33:54 +0000 1132) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1133)         if shell_size < max(d1 + single1.pid_particle_pair[1].radius *
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1134)                            (1.0 + self.SINGLE_SHELL_FACTOR), \
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1135)                                d2 + single2.pid_particle_pair[1].radius * \
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1136)                                (1.0 + self.SINGLE_SHELL_FACTOR)) * 1.3:
049b6bdd (mozo                    2009-01-26 06:54:38 +0000 1137)             if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1138)                 log.debug('%s not formed: singles are better' %
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1139)                       'Pair(%s, %s)' % (single1.pid_particle_pair[0], 
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1140)                                           single2.pid_particle_pair[0]))
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1141)             return None
e3e19d34 (shafi                   2007-11-28 23:28:03 +0000 1142) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1143)         # 3. Ok, Pair makes sense.  Create one.
7790978e (thomie                  2010-03-19 20:27:36 +0900 1144)         shell_size = min(shell_size, max_shell_size)
e3e19d34 (shafi                   2007-11-28 23:28:03 +0000 1145) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1146)         pair = self.create_pair(single1, single2, com, r0, shell_size)
7e3d1fd5 (thomie                  2010-02-23 15:27:07 +0900 1147) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1148)         pair.dt, pair.event_type, pair.reactingsingle = \
7790978e (thomie                  2010-03-19 20:27:36 +0900 1149)             pair.determine_next_event(r0)
c800ba97 (thomie                  2010-02-24 14:43:36 +0900 1150) 
7e3d1fd5 (thomie                  2010-02-23 15:27:07 +0900 1151)         assert pair.dt >= 0
6ecd03ee (shafi                   2009-11-26 09:18:06 +0000 1152) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1153)         self.last_time = self.t
d8b8497a (shafi                   2007-11-26 11:10:22 +0000 1154) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1155)         self.remove_domain(single1)
7790978e (thomie                  2010-03-19 20:27:36 +0900 1156)         self.remove_domain(single2)
117beb87 (shafi                   2007-09-13 00:53:41 +0000 1157) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1158)         self.add_pair_event(pair)
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1159)         # single1 will be removed by the scheduler.
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1160)         self.removeEvent(single2)
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1161) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1162)         assert closest_shell_distance == numpy.inf or shell_size < closest_shell_distance
7790978e (thomie                  2010-03-19 20:27:36 +0900 1163)         assert shell_size >= min_shell_size_with_margin
7790978e (thomie                  2010-03-19 20:27:36 +0900 1164)         assert shell_size <= max_shell_size
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1165) 
049b6bdd (mozo                    2009-01-26 06:54:38 +0000 1166)         if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1167)             log.info('%s, dt=%g, r0=%g, shell=%g,' %
7790978e (thomie                  2010-03-19 20:27:36 +0900 1168)                  (pair, pair.dt, r0, shell_size) + 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1169)                  'closest=%s, closest_shell_distance=%g' %
7790978e (thomie                  2010-03-19 20:27:36 +0900 1170)                  (closest, closest_shell_distance))
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1171) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1172)         assert self.check_obj(pair)
11c410e8 (shafi                   2008-06-09 01:10:56 +0000 1173) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1174)         return pair
84b7b61b (shafi                   2007-04-26 02:05:20 +0000 1175)     
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1176) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1177)     def form_multi(self, single, neighbors, dists):
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1178) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1179)         min_shell = single.pid_particle_pair[1].radius * (1.0 + self.MULTI_SHELL_FACTOR)
31258958 (thomie                  2010-02-05 14:41:35 +0900 1180)         # Multis shells need to be contiguous.
7790978e (thomie                  2010-03-19 20:27:36 +0900 1181)         if dists[0] > min_shell:
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1182)             return None
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1183) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1184)         neighbors = [neighbors[i] for i in (dists <= min_shell).nonzero()[0]]
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1185) 
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1186)         closest = neighbors[0]
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1187) 
1934056a (shafi                   2009-12-18 07:44:39 +0000 1188)         # if the closest to this Single is a Single, create a new Multi
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1189)         if isinstance(closest, Single):
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1190) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1191)             multi = self.create_multi()
7790978e (thomie                  2010-03-19 20:27:36 +0900 1192)             self.add_to_multi(single, multi)
7790978e (thomie                  2010-03-19 20:27:36 +0900 1193)             self.remove_domain(single)
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1194)             for neighbor in neighbors:
7790978e (thomie                  2010-03-19 20:27:36 +0900 1195)                 self.add_to_multi_recursive(neighbor, multi)
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1196) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1197)             multi.initialize(self.t)
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1198)             
7790978e (thomie                  2010-03-19 20:27:36 +0900 1199)             self.add_multi_event(multi)
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1200) 
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1201)             return multi
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1202) 
1934056a (shafi                   2009-12-18 07:44:39 +0000 1203)         # if the closest to this Single is a Multi, reuse the Multi.
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1204)         elif isinstance(closest, Multi):
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1205) 
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1206)             multi = closest
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1207)             if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1208)                 log.info('multi merge %s %s' % (single, multi))
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1209) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1210)             self.add_to_multi(single, multi)
7790978e (thomie                  2010-03-19 20:27:36 +0900 1211)             self.remove_domain(single)
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1212)             for neighbor in neighbors[1:]:
7790978e (thomie                  2010-03-19 20:27:36 +0900 1213)                 self.add_to_multi_recursive(neighbor, multi)
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1214) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1215)             multi.initialize(self.t)
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1216) 
56803120 (Moriyoshi Koizumi       2010-06-23 15:58:42 +0900 1217)             self.update_multi_event(self.t + multi.dt, multi)
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1218) 
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1219)             return multi
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1220) 
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1221) 
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1222)         assert False, 'do not reach here'
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1223) 
ce668d73 (shafi                   2009-12-18 07:42:24 +0000 1224) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1225)     def add_to_multi_recursive(self, obj, multi):
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1226)         if isinstance(obj, Single):
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1227)             if multi.has_particle(obj.pid_particle_pair[0]):  # Already in the Multi.
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1228)                 return
7790978e (thomie                  2010-03-19 20:27:36 +0900 1229)             assert obj.is_reset()
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900 1230)             objpos = obj.shell.shape.position
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1231)             
7790978e (thomie                  2010-03-19 20:27:36 +0900 1232)             self.add_to_multi(obj, multi)
7790978e (thomie                  2010-03-19 20:27:36 +0900 1233)             self.remove_domain(obj)
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1234)             self.removeEvent(obj)
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1235) 
2ca8ed51 (moriyoshi               2009-12-09 04:32:43 +0000 1236)             radius = obj.pid_particle_pair[1].radius *\
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1237)                 (1.0 + self.MULTI_SHELL_FACTOR)
7790978e (thomie                  2010-03-19 20:27:36 +0900 1238)             neighbors = self.get_neighbors_within_radius_no_sort(objpos, radius,
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1239)                                                             ignore=[obj.domain_id])
c9813def (shafi                   2009-12-09 09:58:21 +0000 1240) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1241)             burst = self.burst_non_multis(neighbors)
7790978e (thomie                  2010-03-19 20:27:36 +0900 1242)             neighbor_dists = self.obj_distance_array(objpos, burst)
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1243)             neighbors = [burst[i] for i in 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1244)                          (neighbor_dists <= radius).nonzero()[0]]
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1245) 
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1246)             for obj in neighbors:
7790978e (thomie                  2010-03-19 20:27:36 +0900 1247)                 self.add_to_multi_recursive(obj, multi)
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1248) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1249)         elif isinstance(obj, Multi):
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1250)             for pp in multi.particles:
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1251)                 if obj.has_particle(pp[0]):
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1252)                     if __debug__:
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1253)                         log.debug('%s already added. skipping.' % obj)
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1254)                     break
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1255)             else:
7790978e (thomie                  2010-03-19 20:27:36 +0900 1256)                 self.merge_multis(obj, multi)
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1257)         else:
ff8f1af9 (shafi                   2009-12-10 10:06:48 +0000 1258)             assert False, 'do not reach here.'  # Pairs are burst
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1259) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1260)     def new_spherical_shell(self, domain_id, pos, size):
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1261)         shell_id_shell_pair = (
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1262)             self.shell_id_generator(),
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1263)             SphericalShell(domain_id, Sphere(pos, size)))
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1264)         self.move_shell(shell_id_shell_pair)
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1265)         return shell_id_shell_pair
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1266) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1267)     def add_to_multi(self, single, multi):
049b6bdd (mozo                    2009-01-26 06:54:38 +0000 1268)         if __debug__:
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1269)             log.info('adding %s to %s' % (single.domain_id, multi.domain_id))
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1270)             log.info(single)
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1271)             log.info(multi)
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1272)         sid_shell_pair = self.new_spherical_shell(
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1273)             multi.domain_id,
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1274)             single.pid_particle_pair[1].position,
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1275)             single.pid_particle_pair[1].radius * \
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1276)                 (1.0 + self.MULTI_SHELL_FACTOR))
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1277)         multi.add_shell(sid_shell_pair)
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1278)         multi.add_particle(single.pid_particle_pair)
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1279) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1280)     def merge_multis(self, multi1, multi2):
60be4c82 (moriyoshi               2009-05-01 01:52:28 +0000 1281)         '''
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1282)         merge multi1 into multi2. multi1 will be removed.
60be4c82 (moriyoshi               2009-05-01 01:52:28 +0000 1283)         '''
049b6bdd (mozo                    2009-01-26 06:54:38 +0000 1284)         if __debug__:
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1285)             log.info('merging %s to %s' % (multi1.domain_id, multi2.domain_id))
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1286)             log.info(multi1)
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1287)             log.info(multi2)
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1288) 
c87be851 (Moriyoshi Koizumi       2010-03-16 15:52:04 +0900 1289)             try:
e2dbe6de (Moriyoshi Koizumi       2010-07-05 17:48:26 +0900 1290)                 some_particle_of_multi1 = iter(multi1.particle_container).next()
e2dbe6de (Moriyoshi Koizumi       2010-07-05 17:48:26 +0900 1291)                 assert some_particle_of_multi1[0] not in multi2.particle_container
c87be851 (Moriyoshi Koizumi       2010-03-16 15:52:04 +0900 1292)             except:
c87be851 (Moriyoshi Koizumi       2010-03-16 15:52:04 +0900 1293)                 pass
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1294) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1295)         for sid_shell_pair in multi1.shell_list:
7a5f0206 (Koichi Takahashi        2010-05-13 17:11:10 +0900 1296)             sid_shell_pair[1].did = multi2.domain_id
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1297)             self.move_shell(sid_shell_pair)
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1298)             multi2.add_shell(sid_shell_pair)
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1299) 
7a5f0206 (Koichi Takahashi        2010-05-13 17:11:10 +0900 1300)         for pid_particle_pair in multi1.particles:
7a5f0206 (Koichi Takahashi        2010-05-13 17:11:10 +0900 1301)             multi2.add_particle(pid_particle_pair)
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1302) 
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1303)         del self.domains[multi1.domain_id]
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1304)         self.removeEvent(multi1)
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1305) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1306)     def get_neighbors_within_radius_no_sort(self, pos, radius, ignore=[]):
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1307)         """Get neighbor domains within given radius.
e4c5fde9 (shafi                   2009-12-10 08:03:44 +0000 1308) 
e4c5fde9 (shafi                   2009-12-10 08:03:44 +0000 1309)         ignore: domain ids.
fd242469 (shafi                   2007-07-04 02:42:21 +0000 1310) 
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1311)         Only returns neighbors, not the distances towards their shells. Can 
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1312)         for example be used to try to clear all objects from a certain volume.
3d8c644d (shafi                   2009-12-10 09:44:06 +0000 1313) 
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1314)         """
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1315)         neighbors = []
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1316)         for container in self.containers:
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1317)             result = container.get_neighbors_within_radius(pos, radius)
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1318)             # result = [((shell_id_shell_pair), distance), ]
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1319)             # Since a domain can have more than 1 shell (multis for example), 
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1320)             # and for each shell there is an entry in the shell container, we 
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1321)             # make sure each domain occurs only once in the returned list 
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1322)             # here.
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1323)             neighbors.extend(self.domains[did]
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1324)                              for did in uniq(s[0][1].did for s in result)
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1325)                                      if did not in ignore)
47d8c17e (thomie                  2010-02-26 22:35:42 +0900 1326)         return neighbors
fd242469 (shafi                   2007-07-04 02:42:21 +0000 1327) 
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1328)     def get_intruders(self, position, radius, ignore):
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1329)         intruders = []   # intruders are domains within radius
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1330)         closest_domain = None   # closest domain, excluding intruders.
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1331)         closest_distance = numpy.inf # distance to the shell of the closest.
3e16b161 (shafi                   2008-04-25 04:14:29 +0000 1332) 
8f12364a (moriyoshi               2009-12-09 04:23:51 +0000 1333)         seen = set(ignore)
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1334)         for container in self.containers:
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1335)             neighbors = container.get_neighbors(position)
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1336)             for n in neighbors:
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1337)                 domain_id = n[0][1].did
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1338)                 distance = n[1]
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1339)                 if distance > radius:
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1340)                     if distance < closest_distance:
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1341)                         # This is domain (the first one for this container) 
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1342)                         # that has a shell that is more than radius away from 
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1343)                         # pos.  If it is closer than the closest such one we 
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1344)                         # found so far: store it. Always break out of the 
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1345)                         # inner for loop and check the other containers.
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1346)                         closest_domain = self.domains[domain_id]
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1347)                         closest_distance = distance
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1348)                         break
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1349)                     else:
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1350)                         break
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1351)                 elif domain_id not in seen:
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1352)                     # Since a domain can have more than 1 shell (multis for 
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1353)                     # example), and for each shell there is an entry in the 
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1354)                     # shell container, we make sure each domain occurs only 
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1355)                     # once in the returned list here.
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1356)                     seen.add(domain_id)
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1357)                     intruders.append(self.domains[domain_id])
3e16b161 (shafi                   2008-04-25 04:14:29 +0000 1358) 
04268e43 (thomie                  2010-02-27 00:01:31 +0900 1359)         return intruders, closest_domain, closest_distance
39da29dc (shafi                   2008-04-25 06:02:00 +0000 1360) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1361)     def get_closest_obj(self, pos, ignore=[]):
c9813def (shafi                   2009-12-09 09:58:21 +0000 1362)         '''
c9813def (shafi                   2009-12-09 09:58:21 +0000 1363)         ignore: domain ids.
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1364) 
c9813def (shafi                   2009-12-09 09:58:21 +0000 1365)         '''
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1366)         closest_domain = None
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1367)         closest_distance = numpy.inf
c9813def (shafi                   2009-12-09 09:58:21 +0000 1368) 
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1369)         for container in self.containers:
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1370)             result = container.get_neighbors(pos)
691049b1 (shafi                   2008-04-29 23:19:04 +0000 1371) 
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1372)             for shell_id_shell_pair, distance in result:
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1373)                 domain_id = shell_id_shell_pair[1].did 
691049b1 (shafi                   2008-04-29 23:19:04 +0000 1374) 
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1375)                 if domain_id not in ignore and distance < closest_distance:
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1376)                     domain = self.domains[domain_id]
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1377)                     closest_domain, closest_distance = domain, distance
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1378)                     # Found yet a closer domain. Break out of inner for loop 
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1379)                     # and check other containers.
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1380)                     break   
691049b1 (shafi                   2008-04-29 23:19:04 +0000 1381) 
03ec97f5 (thomie                  2010-02-26 23:24:44 +0900 1382)         return closest_domain, closest_distance
691049b1 (shafi                   2008-04-29 23:19:04 +0000 1383) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1384)     def obj_distance(self, pos, obj):
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900 1385)         return min(self.world.distance(shell.shape, pos) for i, (_, shell) in enumerate(obj.shell_list))
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1386) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1387)     def obj_distance_array(self, pos, objs):
7790978e (thomie                  2010-03-19 20:27:36 +0900 1388)         dists = numpy.array([self.obj_distance(pos, obj) for obj in objs])
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1389)         return dists
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1390)             
75e8dfa0 (shafi                   2007-06-30 04:22:03 +0000 1391) 
df083205 (shafi                   2007-09-19 23:30:15 +0000 1392)     #
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1393)     # statistics reporter
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1394)     #
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1395) 
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1396)     def print_report(self, out=None):
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1397)         report = '''
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1398) t = %g
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1399) steps = %d 
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1400) \tSingle:\t%d\t(escape: %d, reaction: %d)
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1401) \tPair:\t%d\t(escape r: %d, R: %d, reaction pair: %d, single: %d)
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1402) \tMulti:\t%d\t(escape: %d, reaction: %d)
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1403) total reactions = %d
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1404) rejected moves = %d
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1405) '''\
7790978e (thomie                  2010-03-19 20:27:36 +0900 1406)             % (self.t, self.step_counter,
dc99cd59 (shafi                   2009-12-15 06:03:27 +0000 1407)                numpy.array(self.single_steps.values()).sum(),
616a0341 (thomie                  2010-02-22 17:13:13 +0900 1408)                self.single_steps[EventType.SINGLE_ESCAPE],
616a0341 (thomie                  2010-02-22 17:13:13 +0900 1409)                self.single_steps[EventType.SINGLE_REACTION],
dc99cd59 (shafi                   2009-12-15 06:03:27 +0000 1410)                numpy.array(self.pair_steps.values()).sum(),
616a0341 (thomie                  2010-02-22 17:13:13 +0900 1411)                self.pair_steps[EventType.IV_ESCAPE],
616a0341 (thomie                  2010-02-22 17:13:13 +0900 1412)                self.pair_steps[EventType.COM_ESCAPE],
91a21f9f (thomie                  2010-02-25 20:49:52 +0900 1413)                self.pair_steps[EventType.IV_REACTION],
616a0341 (thomie                  2010-02-22 17:13:13 +0900 1414)                self.pair_steps[EventType.SINGLE_REACTION],
f86d6843 (shafi                   2009-12-21 07:21:26 +0000 1415)                self.multi_steps[2], # total multi steps
616a0341 (thomie                  2010-02-22 17:13:13 +0900 1416)                self.multi_steps[EventType.MULTI_ESCAPE],
616a0341 (thomie                  2010-02-22 17:13:13 +0900 1417)                self.multi_steps[EventType.MULTI_REACTION],
7790978e (thomie                  2010-03-19 20:27:36 +0900 1418)                self.reaction_events,
7790978e (thomie                  2010-03-19 20:27:36 +0900 1419)                self.rejected_moves
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1420)                )
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1421) 
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1422)         print >> out, report
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1423) 
0c22362f (shafi                   2009-12-15 05:41:14 +0000 1424)     #
df083205 (shafi                   2007-09-19 23:30:15 +0000 1425)     # consistency checkers
df083205 (shafi                   2007-09-19 23:30:15 +0000 1426)     #
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1427) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1428)     def check_obj(self, obj):
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1429)         obj.check()
c973286b (shafi                   2007-10-17 08:55:55 +0000 1430) 
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000 1431)         for shell_id, shell in obj.shell_list:
7790978e (thomie                  2010-03-19 20:27:36 +0900 1432)             closest, distance = self.get_closest_obj(shell.shape.position,
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1433)                                                    ignore = [obj.domain_id])
737c6f80 (thomie                  2010-03-13 00:19:25 +0900 1434)             if(type(obj) is CylindricalSurfaceSingle or
737c6f80 (thomie                  2010-03-13 00:19:25 +0900 1435)                type(obj) is CylindricalSurfacePair):
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900 1436)                 shell_size = shell.shape.size
bf8ee670 (thomie                  2010-02-26 22:08:33 +0900 1437)             else:
97e031f6 (Moriyoshi Koizumi       2010-03-17 15:53:46 +0900 1438)                 shell_size = shell.shape.radius
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1439) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1440)             assert shell_size <= self.get_user_max_shell_size(),\
11c410e8 (shafi                   2008-06-09 01:10:56 +0000 1441)                 '%s shell size larger than user-set max shell size' % \
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1442)                 str(shell_id)
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1443) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1444)             assert shell_size <= self.get_max_shell_size(),\
11c410e8 (shafi                   2008-06-09 01:10:56 +0000 1445)                 '%s shell size larger than simulator cell size / 2' % \
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1446)                 str(shell_id)
f87104f5 (shafi                   2008-06-08 04:12:33 +0000 1447) 
bf8ee670 (thomie                  2010-02-26 22:08:33 +0900 1448)             assert distance - shell_size >= 0.0,\
11c410e8 (shafi                   2008-06-09 01:10:56 +0000 1449)                 '%s overlaps with %s. (shell: %g, dist: %g, diff: %g.' \
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1450)                 % (str(obj), str(closest), shell_size, distance,\
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1451)                        distance - shell_size)
faba9096 (shafi                   2008-04-22 05:17:09 +0000 1452) 
11c410e8 (shafi                   2008-06-09 01:10:56 +0000 1453)         return True
75e8dfa0 (shafi                   2007-06-30 04:22:03 +0000 1454) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1455)     def check_obj_for_all(self):
bed38b09 (Moriyoshi Koizumi       2010-07-05 16:21:47 +0900 1456)         for id, event in self.scheduler:
bed38b09 (Moriyoshi Koizumi       2010-07-05 16:21:47 +0900 1457)             self.check_obj(event.data)
5a5dacbb (shafi                   2007-09-24 23:14:31 +0000 1458) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1459)     def check_event_stoichiometry(self):
7790978e (thomie                  2010-03-19 20:27:36 +0900 1460)         event_population = 0
bed38b09 (Moriyoshi Koizumi       2010-07-05 16:21:47 +0900 1461)         for id, event in self.scheduler:
bed38b09 (Moriyoshi Koizumi       2010-07-05 16:21:47 +0900 1462)             event_population += event.data.multiplicity
5a5dacbb (shafi                   2007-09-24 23:14:31 +0000 1463) 
84afe813 (Moriyoshi Koizumi       2010-04-26 20:45:20 +0900 1464)         if self.world.num_particles != event_population:
7790978e (thomie                  2010-03-19 20:27:36 +0900 1465)             raise RuntimeError, 'population %d != event_population %d' %\
7790978e (thomie                  2010-03-19 20:27:36 +0900 1466)                   (population, event_population)
dd6ae27b (shafi                   2007-11-27 00:09:18 +0000 1467) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1468)     def check_shell_matrix(self):
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1469)         did_map = {}
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1470)         shell_map = {}
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900 1471)         for container in self.containers:
acd00e2c (Moriyoshi Koizumi       2010-03-11 14:21:35 +0900 1472)             if self.world.world_size != container.world_size:
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900 1473)                 raise RuntimeError,\
acd00e2c (Moriyoshi Koizumi       2010-03-11 14:21:35 +0900 1474)                     'self.world.world_size != container.world_size'
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1475)             for shell_id, shell in container:
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1476)                 did_map.setdefault(shell.did, []).append(shell_id)
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1477)                 shell_map[shell_id] = shell
d6f015d5 (shafi                   2007-12-12 17:09:35 +0000 1478) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1479)         shell_population = 0
bed38b09 (Moriyoshi Koizumi       2010-07-05 16:21:47 +0900 1480)         for id, event in self.scheduler:
bed38b09 (Moriyoshi Koizumi       2010-07-05 16:21:47 +0900 1481)             shell_population += event.data.num_shells
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1482)             shell_ids = did_map[event.data.domain_id]
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1483)             if len(shell_ids) != event.data.num_shells:
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1484)                 diff = set(sid for sid, _ in event.data.shell_list).difference(shell_ids)
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1485)                 for sid in diff:
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1486)                     print shell_map.get(sid, None)
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1487) 
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1488)                 raise RuntimeError("number of shells are inconsistent (%d != %d; %s) - %s" % (len(shell_ids), event.data.num_shells, event.data.domain_id, diff))
d08fc8eb (Moriyoshi Koizumi       2010-07-05 20:03:51 +0900 1489) 
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900 1490)         matrix_population = sum(len(container) for container in self.containers)
7790978e (thomie                  2010-03-19 20:27:36 +0900 1491)         if shell_population != matrix_population:
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900 1492)             raise RuntimeError(
d6dc31ff (thomie                  2010-02-26 20:51:09 +0900 1493)                 'num shells (%d) != matrix population (%d)' % 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1494)                 (shell_population, matrix_population))
fe3d6de6 (shafi                   2008-03-18 03:09:45 +0000 1495) 
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1496)     def check_domains(self):
c7ade95b (Moriyoshi Koizumi       2010-07-05 18:41:31 +0900 1497)         event_ids = set(domain.event_id for domain in self.domains.itervalues())
bed38b09 (Moriyoshi Koizumi       2010-07-05 16:21:47 +0900 1498)         for id, event in self.scheduler:
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1499)             if id not in event_ids:
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1500)                 raise RuntimeError,\
bed38b09 (Moriyoshi Koizumi       2010-07-05 16:21:47 +0900 1501)                     '%s in EventScheduler not in self.domains' % event.data
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1502)             event_ids.remove(id)
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1503) 
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1504)         # self.domains always include a None  --> this can change in future
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1505)         if event_ids:
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1506)             raise RuntimeError,\
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1507)                 'following domains in self.domains not in Event Scheduler: %s' \
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1508)                 % str(tuple(event_ids))
380a7581 (shafi                   2008-03-21 03:04:42 +0000 1509) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1510)     def check_pair_pos(self, pair, pos1, pos2, com, radius):
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000 1511)         particle1 = pair.single1.pid_particle_pair[1]
b159dcf1 (moriyoshi               2009-12-04 13:20:11 +0000 1512)         particle2 = pair.single2.pid_particle_pair[1]
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1513) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1514)         old_com = com
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1515)         
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1516)         # debug: check if the new positions are valid:
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900 1517)         new_distance = self.world.distance(pos1, pos2)
7790978e (thomie                  2010-03-19 20:27:36 +0900 1518)         particle_radius12 = particle1.radius + particle2.radius
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1519) 
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1520)         # check 1: particles don't overlap.
7790978e (thomie                  2010-03-19 20:27:36 +0900 1521)         if new_distance <= particle_radius12:
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1522)             if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1523)                 log.info('rejected move: radii %g, particle distance %g',
7790978e (thomie                  2010-03-19 20:27:36 +0900 1524)                          (particle1.radius + particle2.radius, new_distance))
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1525)             if __debug__:
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1526)                 log.debug('DEBUG: pair.dt %g, pos1 %s, pos2 %s' %
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1527)                           (pair.dt, str(pos1), str(pos2)))
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1528)             raise RuntimeError, 'New particles overlap'
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1529) 
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1530)         # check 2: particles within mobility radius.
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900 1531)         d1 = self.world.distance(old_com, pos1) + particle1.radius
8617a12c (Moriyoshi Koizumi       2010-05-20 15:14:55 +0900 1532)         d2 = self.world.distance(old_com, pos2) + particle2.radius
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1533)         if d1 > radius or d2 > radius:
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1534)             raise RuntimeError, \
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1535)                 'New particle(s) out of protective sphere. %s' % \
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1536)                 'radius = %g, d1 = %g, d2 = %g ' % (radius, d1, d2)
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1537)                 
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1538)         
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1539) 
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1540)         return True
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1541) 
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1542) 
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1543) 
7bcccd6f (shafi                   2009-11-26 12:36:17 +0000 1544) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1545)     def check(self):
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1546)         ParticleSimulatorBase.check(self)
380a7581 (shafi                   2008-03-21 03:04:42 +0000 1547) 
85d8e4dd (shafi                   2008-04-24 02:09:01 +0000 1548)         assert self.scheduler.check()
85d8e4dd (shafi                   2008-04-24 02:09:01 +0000 1549) 
beec6e3b (shafi                   2007-04-04 03:34:17 +0000 1550)         assert self.t >= 0.0
beec6e3b (shafi                   2007-04-04 03:34:17 +0000 1551)         assert self.dt >= 0.0
d6f015d5 (shafi                   2007-12-12 17:09:35 +0000 1552) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1553)         self.check_shell_matrix()
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1554)         self.check_domains()
7790978e (thomie                  2010-03-19 20:27:36 +0900 1555)         self.check_event_stoichiometry()
beec6e3b (shafi                   2007-04-04 03:34:17 +0000 1556)         
7790978e (thomie                  2010-03-19 20:27:36 +0900 1557)         self.check_obj_for_all()
ec1b4735 (shafi                   2007-04-26 04:20:06 +0000 1558) 
df083205 (shafi                   2007-09-19 23:30:15 +0000 1559)     #
df083205 (shafi                   2007-09-19 23:30:15 +0000 1560)     # methods for debugging.
df083205 (shafi                   2007-09-19 23:30:15 +0000 1561)     #
08ce7573 (shafi                   2007-07-04 10:08:21 +0000 1562) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1563)     def dump_scheduler(self):
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1564)         for id, event in self.scheduler:
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1565)             print id, event
08ce7573 (shafi                   2007-07-04 10:08:21 +0000 1566) 
a034b53b (thomie                  2010-03-19 20:09:23 +0900 1567)     def dump(self):
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1568)         for id, event in self.scheduler:
bed38b09 (Moriyoshi Koizumi       2010-07-05 16:21:47 +0900 1569)             print id, event, event.data
08ce7573 (shafi                   2007-07-04 10:08:21 +0000 1570) 
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1571)     def count_domains(self):
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1572)         '''
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1573)         Returns a tuple (# Singles, # Pairs, # Multis).
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1574)         '''
08ce7573 (shafi                   2007-07-04 10:08:21 +0000 1575) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1576)         num_singles = 0
7790978e (thomie                  2010-03-19 20:27:36 +0900 1577)         num_pairs = 0
7790978e (thomie                  2010-03-19 20:27:36 +0900 1578)         num_multis = 0
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1579)         for d in self.domains.itervalues():
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1580)             if isinstance(d, Single):
7790978e (thomie                  2010-03-19 20:27:36 +0900 1581)                 num_singles += 1
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1582)             elif isinstance(d, Pair):
7790978e (thomie                  2010-03-19 20:27:36 +0900 1583)                 num_pairs += 1
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1584)             elif isinstance(d, Multi):
7790978e (thomie                  2010-03-19 20:27:36 +0900 1585)                 num_multis += 1
4a60a2b2 (shafi                   2009-12-09 07:33:53 +0000 1586)             else:
274d79b2 (Koichi Takahashi        2010-05-28 16:08:48 +0900 1587)                 raise RuntimeError, 'DO NOT GET HERE'
08ce7573 (shafi                   2007-07-04 10:08:21 +0000 1588) 
7790978e (thomie                  2010-03-19 20:27:36 +0900 1589)         return (num_singles, num_pairs, num_multis)
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1590) 
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1591)     dispatch = [
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1592)         (Single, fire_single),
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1593)         (Pair, fire_pair),
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1594)         (Multi, fire_multi)
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1595)         ]
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1596) 
8cac304a (Moriyoshi Koizumi       2010-07-05 16:07:44 +0900 1597) 
