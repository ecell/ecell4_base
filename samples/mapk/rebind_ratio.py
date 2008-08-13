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
    KpPartner = {}
    PartnerKp = {}

    for r in reactions:
        print r.t, r.reactants, r.products

        # unbinding
        if len( r.reactants ) == 1 and len( r.products ) == 2:
            for i, p in enumerate( r.products ):
                if p.speciesName == 'Kp':
                    Kp = p
                    peer = r.products[ 1 - i ]
                    if peer.speciesName == 'KKi':
                        print 'unbinding', Kp, peer
                        KpCreated[Kp] = r.t
                        KpPartner[Kp] = peer
                        PartnerKp[peer] = Kp

        # binding
        elif len( r.reactants ) == 2 and len( r.products ) == 1:
            for i, p in enumerate( r.reactants ):
                if p.speciesName == 'Kp':
                    Kp = p
                    peer = r.reactants[ 1 - i ]
                    if peer.speciesName == 'KK':
                        print 'binding', Kp, peer
                        t_u = KpCreated[Kp]
                        t_rebinding = r.t - t_u

                        partner = KpPartner[Kp]
                        if peer == partner:
                            print t_rebinding, '\trebinding'
                        else:
                            print t_rebinding, '\tdiffusion'

                        del KpCreated[Kp]
                        del KpPartner[Kp]
                        del PartnerKp[partner]


        # monomolecular
        elif len( r.reactants ) == 1 and len( r.products ) == 1:
            for p in r.reactants:
                if PartnerKp.has_key( p ):
                    Kp = PartnerKp[p]
                    newpartner = r.products[0]
                    KpPartner[Kp] = newpartner
                    del PartnerKp[p]
                    PartnerKp[newpartner] = Kp
                    #print 'tracked', p, ' -> ', newpartner
                                  
        


if __name__ == '__main__':

    import sys

    for file in sys.argv[1:]:
        file = open( sys.argv[1] )
        reactions = load_reactions( file )
        print >> sys.stderr, 'num reactions: ', len( reactions )
        rebind_ratio( reactions )

