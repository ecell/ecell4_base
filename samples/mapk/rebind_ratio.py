#!/usr/bin/env python

import sys

class Particle:
    def __init__( self, particle ):
        self.speciesName = particle[0]
        self.id = particle[1]

    def __eq__( self, other ):
        return self.speciesName == other.speciesName and self.id == other.id



    def __str__( self ):
        return "( '%s', %d )" % ( self.speciesName, self.id )

    def __repr__( self ):
        return self.__str__()

    def __hash__( self ):
        return hash( self.speciesName ) ^ self.id


class ReactionEvent:
    def __init__( self, t, reactants, products ):
        self.t = t
        self.reactants = reactants
        self.products = products


def load_reactions( file ):

    reactions = []

    for line in file.readlines():
        line = eval( line )
        t = line[0]
        reactants = [ Particle( p ) for p in line[1] ]
        products = [ Particle( p ) for p in line[2] ]
        reactions.append( ReactionEvent( t, reactants, products ) )

    return reactions

def rebind_ratio( reactions ):

    KpCreated = {}
    KpKK = {}    

    KpCurrentForm = {}
    CurrentFormKp = {}
    KKCurrentForm = {}
    CurrentFormKK = {}

    counter = 0
    lasttime = 0

    for r in reactions:
        #print r.t, r.reactants, r.products

        lasttime = r.t

        # unbinding
        if len( r.reactants ) == 1 and len( r.products ) == 2:

            # K_KK -> Kp + KKi
            # Kp first site phosphorylated.
            if r.reactants[0].speciesName == 'K_KK':

                for i, p in enumerate( r.products ):
                    if p.speciesName == 'Kp':
                        Kp = p
                        KKi = r.products[ 1 - i ]
                        if KKi.speciesName == 'KKi':

                            #print r.reactants[0], '->', Kp,  KKi
                            KpCreated[Kp] = r.t
                            
                            KpKK[Kp] = KKi

                            KpCurrentForm[Kp] = Kp
                            CurrentFormKp[Kp] = Kp
                            KKCurrentForm[KKi] = KKi
                            CurrentFormKK[KKi] = KKi

            elif r.reactants[0].speciesName == 'Kp_KK':

                Kp_KK = r.reactants[0]

                # Kp_KK -> Kp + KK
                for i, p in enumerate( r.products ):
                    if p.speciesName == 'Kp':
                        Kp = p
                        KK = r.products[ 1 - i ]
                        if KK.speciesName == 'KK':

                            #print r.reactants[0], '->', Kp, KK

                            if CurrentFormKp.has_key(Kp_KK):
                                originalKp = CurrentFormKp[Kp_KK]
                                KpCurrentForm[originalKp] = Kp
                                del CurrentFormKp[Kp_KK]
                                CurrentFormKp[Kp] = originalKp

                            if CurrentFormKK.has_key(Kp_KK):
                                originalKK = CurrentFormKK[Kp_KK]
                                KKCurrentForm[originalKK] = KK
                                del CurrentFormKK[Kp_KK]
                                CurrentFormKK[KK] = originalKK

                            break


                # Kp_KK -> Kpp + KKi
                for i, p in enumerate( r.products ):
                    if p.speciesName == 'Kpp':
                        Kpp = p
                        KKi = r.products[ 1 - i ]
                        if KKi.speciesName == 'KKi':

                            #print r.reactants[0], '->', Kpp, KKi

                            originalKp=None
                            originalKK=None

                            if CurrentFormKK.has_key(Kp_KK):
                                originalKK = CurrentFormKK[Kp_KK]
                                del KKCurrentForm[originalKK]
                                del CurrentFormKK[Kp_KK]

                            if CurrentFormKp.has_key( Kp_KK ):

                                originalKp = CurrentFormKp[Kp_KK]
                                del KpCurrentForm[originalKp]
                                del CurrentFormKp[Kp_KK]

                                #print originalKp

#                                 t_create = KpCreated[originalKp]
#                                 t = r.t - t_create
#                                 partner = KpKK[originalKp]


#                                 if originalKK is not None and originalKK == partner:
#                                     outfile.write( '%.18g\trebinding\n' % t )
#                                 else:
#                                     outfile.write( '%.18g\tdiffusion\n' % t )



                            


        # binding
        elif len( r.reactants ) == 2 and len( r.products ) == 1:


            # Kp + KK -> Kp_KK
            for i, p in enumerate( r.reactants ):
                if p.speciesName == 'Kp':
                    Kp = p
                    KK = r.reactants[ 1 - i ]
                    if KK.speciesName == 'KK':
                        Kp_KK = r.products[0]
                        assert Kp_KK.speciesName == 'Kp_KK'

                        #print Kp, KK, '->', Kp_KK
                        
                        if not CurrentFormKp.has_key( Kp ):
                            break

                        originalKp = CurrentFormKp[Kp]
                        KpCurrentForm[originalKp] = Kp_KK
                        del CurrentFormKp[Kp]
                        CurrentFormKp[Kp_KK] = originalKp

                        if CurrentFormKK.has_key(KK):
                            originalKK = CurrentFormKK[KK]
                            KKCurrentForm[originalKK] = Kp_KK
                            del CurrentFormKK[KK]
                            CurrentFormKK[Kp_KK] = originalKK
                        else:
                            originalKK = None

                        if KpCreated.has_key(originalKp):
#                             pass
                            t_create = KpCreated[originalKp]
                            t = r.t - t_create
                            partner = KpKK[originalKp]
                            if originalKK is not None and originalKK == partner:
                                outfile.write( '%.18g\trebinding\n' % t )
                            else:
                                outfile.write( '%.18g\tdiffusion\n' % t )
                            counter += 1
                            del KpCreated[originalKp]
                            del KpKK[originalKp]

                        
                        break 


#             # Kp + P -> Kp_P
#             for i, p in enumerate( r.reactants ):
#                 if p.speciesName == 'Kp':
#                     Kp = p
#                     P = r.reactants[ 1 - i ]
#                     if P.speciesName == 'P':
#                         Kp_P = r.products[0]
#                         assert Kp_P.speciesName == 'Kp_P'

#                         if not CurrentFormKp.has_key( Kp ):
#                             break

#                         originalKp = CurrentFormKp[Kp]
#                         KpCurrentForm[originalKp] = Kp_P
#                         del CurrentFormKp[Kp]
#                         CurrentFormKp[Kp_P] = originalKp

#                         if KpCreated.has_key(originalKp):
# #                             pass
#                             del KpCreated[originalKp]
#                             del KpKK[originalKp]

                        
#                         break 



        # monomolecular
        elif len( r.reactants ) == 1 and len( r.products ) == 1:
            if CurrentFormKK.has_key( r.reactants[0] ):
                originalform = CurrentFormKK[ r.reactants[0] ]
                KKCurrentForm[ originalform ] = r.products[0]
                del CurrentFormKK[ r.reactants[0] ]
                CurrentFormKK[ r.products[0] ] = originalform
                #print 'transition', r.reactants[0], '->', r.products[0]

    import numpy
    nonreactions = [ t for t in KpCreated.values() if lasttime-t > 10]
    print 'reactions: ', counter, 'non-reactions', len(nonreactions)
    for t in nonreactions:
        #if t > 60:
        outfile.write( '%.18g\tno-reaction\n' % numpy.inf )
    
                                  
        


if __name__ == '__main__':

    import sys
    import os
    import glob

    for pattern in sys.argv[1:]:
    
        globpattern = pattern.replace('ALL','*')

        l = os.path.basename( os.path.splitext( pattern )[0] )

        outfilename = l + '.rebind'
        outfile = open( outfilename, 'w' )
        print >> sys.stderr, 'pattern ', l, '\noutfile', outfilename

        filelist = glob.glob( globpattern )


        for file in filelist:
            reactions = load_reactions( open( file ) )
            print >> sys.stderr, 'num reactions: ', len( reactions )
            rebind_ratio( reactions )

        outfile.close()
