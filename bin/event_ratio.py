#!/usr/bin/env python

import numpy

import sys
import glob


filepattern = open( sys.argv[1] )

filepattern.replace( 'ALL', '*' )
filelist = glob.glob( filepattern )
print filelist


events = {}

for line in file.readlines():
    line = line.split()
    eventType = line[1]
    t = float(line[0])
    if events.has_key( eventType ):
        events[eventType].append(t)
    else:
        events[eventType] = [t]

total = numpy.sum( [len(v) for v in events.values() ] )

print 'total', total
for e in events.keys():
    print e, 'ratio', float(len(events[e]))/total, 't_mean', numpy.mean(events[e])
