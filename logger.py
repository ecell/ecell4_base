
import os
import re
#import logging
import numpy

import logging

INF = numpy.inf

log = logging.getLogger('ecell')


class Logger:

    def __init__( self, sim, logname = 'log', directory = 'data',
                  comment='' ):

        self.sim = sim


        self.logname = logname

        self.fileCounter = 0

        self.directory = directory
        try:
            os.mkdir( directory )
        except:
            pass

        self.particleOutInterval = INF

        self.lastTime = 0.0
        self.nextTime = INF

        self.particleOutPattern = re.compile( '' )
        self.prepareTimecourseFile( comment )
        self.writeTimecourse()


    def setInterval( self, interval ):
        self.interval = interval

    def setParticleOutPattern( self, pattern ):
        self.particleOutPattern = re.compile( pattern )

    def getParticleOutPattern( self ):
        return self.particleOutPattern.pattern

    def setParticleOutInterval( self, interval ):
        self.particleOutInterval = interval
        self.lastTime = self.sim.t
        self.nextTime = self.lastTime + self.particleOutInterval

    def prepareTimecourseFile( self, comment ):

        self.timecourseFilename = self.logname + '_tc' + '.dat'
        self.timecourseFile = open( self.directory + os.sep +\
                                    self.timecourseFilename, 'w' )
        self.writeTimecourseComment( comment )

        speciesNameList = '\'' +\
            "\', \'".join( str(i) for i in self.sim.speciesList.keys()  ) + '\''
        columns = '[ \'t\', ' + speciesNameList + ']'
        self.writeTimecourseComment( '@ columns= ' + columns )


    def writeTimecourseComment( self, s ):
        self.timecourseFile.write( '#' + s + '\n' )

    def writeTimecourse( self ):
        data = [
            ]

        self.timecourseFile.write( '%g' % self.sim.t + '\t' )
        self.timecourseFile.write( '\t'.join(
            str(len(self.sim.getParticlePool(i.type.id))) \
            for i in self.sim.getSpecies()) + '\n' )
        self.timecourseFile.flush()

    def writeParticles( self ):
        filename = self.logname + '_' + \
            str( self.fileCounter ).zfill(4) + '.dat'

        file = open( self.directory + os.sep + filename, 'w' )

        file.write( '#@ name = \'%s\'\n' % str( self.logname ) )
        file.write( '#@ count = %d\n' % int( self.fileCounter ) )
        file.write( '#@ t = %s\n' % '%g' % self.sim.t )
        file.write( '#@ worldSize = %f\n' % float( self.sim.getWorldSize() ) )
        file.write( '#--------\n' )

        for sid in self.sim.speciesList.keys():
            pid_list = self.sim.particlePool[ sid ]
            for i in pid_list:
                particle = self.sim.particleMatrix[i]
                species = self.sim.speciesList[ sid ]
                file.write( '%s\t%20.14g %20.14g %20.14g %.15g\n' % 
                            ( species.id, particle.position[0], particle.position[1], particle.position[2], species.radius ) )

            file.write( '#\n' )

        file.close()

        self.fileCounter += 1


    def log( self ):

        self.logTimeCourse()
        self.logParticles()

    def logTimeCourse( self ):

        if self.sim.lastReaction:
            self.writeTimecourse()


    def logParticles( self ):
        sim = self.sim
        if self.nextTime <= sim.t + sim.dt:
            #if __debug__: log.info( 'log %g' % self.nextTime )

            sim.stop( self.nextTime )
            self.writeParticles()

            self.nextTime += self.particleOutInterval


        

        
