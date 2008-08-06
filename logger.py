
import os
import re
#import logging
import numpy

INF = numpy.inf


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
            "\', \'".join( self.sim.speciesList.keys()  ) + '\''
        columns = '[ \'t\', ' + speciesNameList + ']'
        self.writeTimecourseComment( '@ columns= ' + columns )


    def writeTimecourseComment( self, s ):
        self.timecourseFile.write( '#' + s + '\n' )

    def writeTimecourse( self ):

        data = [ str( i.pool.size )\
                 for i in self.sim.speciesList.values() ]
            
        self.timecourseFile.write( '%g' % self.sim.t + '\t' )
        self.timecourseFile.write( '\t'.join( data ) + '\n' )
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

        for speciesName in self.sim.speciesList.keys():
            species = self.sim.speciesList[ speciesName ]
            for i in species.pool.positions:
                file.write( '%s\t%20.14g %20.14g %20.14g %.15g\n' % 
                            ( speciesName, i[0], i[1], i[2], species.radius ) )

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
            #log.info( 'log %g' % self.nextTime )

            sim.stop( self.nextTime )
            self.writeParticles()

            self.nextTime += self.particleOutInterval


        

        
