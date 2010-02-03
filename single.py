class Single( object ):
    def __init__( self, domain_id, pid_particle_pair, shell_id_shell_pair, reactiontypes ):
        self.multiplicity = 1

        self.pid_particle_pair = pid_particle_pair
        self.reactiontypes = reactiontypes

        self.k_tot = 0

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.shell_list = [shell_id_shell_pair, ]

        self.eventID = None

        self.domain_id = domain_id

        self.updatek_tot()

    def getD( self ):
        return self.pid_particle_pair[1].D
    D = property( getD )

    def getMinRadius(self):
        return self.pid_particle_pair[1].radius
    minRadius = property(getMinRadius)

    def getShell(self):
        return self.shell_list[0]

    def setShell(self, value):
        self.shell_list[0] = value

    shell = property(getShell, setShell)

    def initialize( self, t ):
        '''
        Initialize this Single.

        The radius (shell size) is shrunken to the radius of the particle
        it represents.   
        self.lastTime is reset to the current time, and self.dt
        is set to zero.
        '''
        self.reset()
        self.lastTime = t

    def reset( self ):
        '''
        Reset the Single.

        Radius (shell size) is shrunken to the actual radius of the particle.
        self.dt is reset to 0.0.  Do not forget to reschedule this Single
        after calling this method.
        '''
        self.dt = 0.0
        self.eventType = EventType.ESCAPE

    def isReset( self ):
        return self.dt == 0.0 and self.eventType == EventType.ESCAPE
        
    def drawR( self, dt, shellSize):
        assert dt >= 0
        a = shellSize - self.pid_particle_pair[1].radius
        rnd = myrandom.uniform()
        gf = FirstPassageGreensFunction( self.pid_particle_pair[1].D )
        gf.seta(a)
        try:
            r = gf.drawR(rnd , dt)
        except Exception, e:
            raise Exception, 'gf.drawR failed; %s; rnd=%g, t=%g, %s' %\
                (str(e), rnd, dt, gf.dump())
        return r

    def drawReactionTime( self ):
        if self.k_tot == 0:
            return numpy.inf
        if self.k_tot == numpy.inf:
            return 0.0
        rnd = myrandom.uniform()
        dt = ( 1.0 / self.k_tot ) * math.log( 1.0 / rnd )
        return dt

    def drawEscapeTime(self, a):
        D = self.pid_particle_pair[1].D

        if D == 0:
            return numpy.inf

        gf = FirstPassageGreensFunction(D)
        gf.seta(a)

        rnd = myrandom.uniform()
        try:
            return gf.drawTime( rnd )
        except Exception, e:
            raise Exception, 'gf.drawTime() failed; %s; rnd=%g, %s' %\
                ( str( e ), rnd, gf.dump() )


    def updatek_tot( self ):
        self.k_tot = 0

        if not self.reactiontypes:
            return

        for rt in self.reactiontypes:
            self.k_tot += rt.k


    def drawReactionRule( self ):
        k_array = [ rt.k for rt in self.reactiontypes ]
        k_array = numpy.add.accumulate( k_array )
        k_max = k_array[-1]

        rnd = myrandom.uniform()
        i = numpy.searchsorted( k_array, rnd * k_max )

        return self.reactiontypes[i]


    def check( self ):
        pass

    def __repr__( self ):
        return 'Single[%s: %s: eventID=%s]' % ( self.domain_id, self.pid_particle_pair[0], self.eventID )


# def calculatePairCoM( pos1, pos2, D1, D2, worldSize ):
#     '''
#     Calculate and return the center-of-mass of a Pair.
#     '''
#     #pos2t = cyclicTranspose( pos2, pos1, worldSize )
#     pos2t = cyclic_transpose( pos2, pos1, worldSize )
#     return ( ( D2 * pos1 + D1 * pos2t ) / ( D1 + D2 ) ) % worldSize



