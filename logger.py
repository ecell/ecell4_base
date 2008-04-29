
import os
import string
import logging


class Logger:

    def __init__( self, simulator, logname = 'log', directory = 'data',
                  comment='' ):

        self.simulator = simulator

        self.interval = 0.001

        self.lastTime = 0.0
        self.nextTime = 0.0

        self.logname = logname

        self.fileCounter = 0

        self.directory = directory
        try:
            os.mkdir( directory )
        except:
            pass

        self.particleOutList = []
        self.prepareTimecourseFile( comment )
        self.writeTimecourse()


    def setInterval( self, interval ):
        self.interval = interval

    def setParticleOutput( self, outlist ):
        self.particleOutList = outlist

    def prepareTimecourseFile( self, comment ):

        self.timecourseFilename = self.logname + '.dat'
        self.timecourseFile = open( self.directory + os.sep +\
                                    self.timecourseFilename, 'w' )
        self.writeTimecourseComment( comment )

        speciesNameList = '\'' +\
            string.join( self.simulator.speciesList.keys(), '\', \'' ) + '\''
        columns = '[ \'t\', ' + speciesNameList + ']'
        self.writeTimecourseComment( '@ columns= ' + columns )


    def writeTimecourseComment( self, s ):
        self.timecourseFile.write( '#' + s + '\n' )

    def writeTimecourse( self ):

        data = [ str( i.pool.size )\
                 for i in self.simulator.speciesList.values() ]
            
        self.timecourseFile.write( str( self.simulator.t ) + '\t' )
        self.timecourseFile.write( string.join( data, '\t' ) + '\n' )
        self.timecourseFile.flush()

    def writeParticles( self ):

        for speciesName in self.particleOutList:
            species = self.simulator.speciesList[ speciesName ]
            positions = species.pool.positions

            filename = speciesName + '_' + '%04d' % self.fileCounter + '.dat'
            #print filename
            file = open( self.directory + os.sep + filename, 'w' )
            file.write( '# name: %s\n' % speciesName )
            file.write( '# radius: %f\n' % species.radius )
            file.write( '# count: %d\n' % species.pool.size )
            file.write( '# t: %f\n' % self.simulator.t )
            file.write( '#--------\n' )

            for i in species.pool.positions:
                file.write( '%.15g\t%.15g\t%.15g\n' % ( i[0],i[1],i[2] ) )

            file.close()

        self.fileCounter += 1


    def log( self ):
        self.logTimeCourse()

        if self.particleOutList:
            self.logParticles()

    def logTimeCourse( self ):

        if self.simulator.populationChanged:
            self.writeTimecourse()

    def logParticles( self ):
        sim = self.simulator
        currentTime = sim.t
        nextTime = sim.t + sim.dt

        if self.nextTime <= nextTime:
            logging.info( 'log %g' % self.nextTime )

            sim.stop( self.nextTime )
            self.writeParticles()

            self.nextTime += self.interval


        

        
