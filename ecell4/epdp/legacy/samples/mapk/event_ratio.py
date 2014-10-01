#!/usr/bin/env python

import numpy

import sys
import glob


file = open(sys.argv[1])

# filepattern.replace('ALL', '*')
# filelist = glob.glob(filepattern)
# print filelist


events = {}


for line in file.readlines():
    line = line.split()
    event_type = line[1]
    t = float(line[0])
    if events.has_key(event_type):
        events[event_type].append(t)
    else:
        events[event_type] = [t]

total = numpy.sum([len(v) for v in events.values()])

print 'total', total
for e in events.keys():
    print e, 'ratio', float(len(events[e]))/total, 't_mean', numpy.mean(events[e])
