from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.species import generate_reactions


@species_attributes
def attributegen():
    egf(r) | 1
    Grb2(SH2, SH3) | 2
    Shc(PTB, Y317=Y) | 3
    Sos(dom) | 4
    egfr(l, r, Y1068=Y, Y1148=Y) | 5
    Grb2(SH2, SH3^1).Sos(dom^1) | 6

@reaction_rules
def rulegen():
    # Ligand-receptor binding (ligand-monomer)
    egfr(l, r) + egf(r) == egfr(l^1, r).egf(r^1) | (1, 2)

    # Note changed multiplicity
    # Receptor-aggregation
    egfr(l^_, r) + egfr(l^_, r) == egfr(l^_,r^3).egfr(l^_,r^3) | (3, 4)
    # egfr(l^_1, r) + egfr(l^_2, r) == egfr(l^_1,r^3).egfr(l^_2,r^3) | (kp2, km2) #XXX: this should work, but not now

    # Transphosphorylation of egfr by RTK
    egfr(r^_, Y1068=Y) > egfr(r^_, Y1068=pY) | 5
    egfr(r^_, Y1148=Y) > egfr(r^_, Y1148=pY) | 6

    # Dephosphorylation
    egfr(Y1068=pY) > egfr(Y1068=Y) | 7
    egfr(Y1148=pY) > egfr(Y1148=Y) | 8

    # Shc transphosphorylation
    egfr(r^_, Y1148=pY^1).Shc(PTB^1,Y317=Y) > egfr(r^_,Y1148=pY^1).Shc(PTB^1,Y317=pY) | 9
    Shc(PTB^_,Y317=pY) > Shc(PTB^_,Y317=Y) | 10

    # Y1068 activity
    egfr(Y1068=pY) + Grb2(SH2,SH3) == egfr(Y1068=pY^1).Grb2(SH2^1,SH3) | (11, 12)
    egfr(Y1068=pY) + Grb2(SH2,SH3^_) == egfr(Y1068=pY^1).Grb2(SH2^1,SH3^_) | (13, 14)
    egfr(Y1068=pY^1).Grb2(SH2^1,SH3) + Sos(dom) == egfr(Y1068=pY^1).Grb2(SH2^1,SH3^2).Sos(dom^2) | (15, 16)

    # Y1148 activity
    egfr(Y1148=pY) + Shc(PTB,Y317=Y) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=Y) | (17, 18)
    egfr(Y1148=pY) + Shc(PTB,Y317=pY) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY) | (19, 20)
    egfr(Y1148=pY) + Shc(PTB,Y317=pY^1).Grb2(SH2^1,SH3) == egfr(Y1148=pY^2).Shc(PTB^2,Y317=pY^1).Grb2(SH2^1,SH3) | (21, 22)
    egfr(Y1148=pY) + Shc(PTB,Y317=pY^1).Grb2(SH2^1,SH3^3).Sos(dom^3) == egfr(Y1148=pY^2).Shc(PTB^2,Y317=pY^1).Grb2(SH2^1,SH3^3).Sos(dom^3) | (23, 24)

    egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY) + Grb2(SH2,SH3) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY^2).Grb2(SH2^2,SH3) | (25, 26)

    egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY) + Grb2(SH2,SH3^3).Sos(dom^3) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY^2).Grb2(SH2^2,SH3^3).Sos(dom^3) | (27, 28)

    Shc(PTB^_,Y317=pY^2).Grb2(SH2^2,SH3) + Sos(dom) == Shc(PTB^_,Y317=pY^2).Grb2(SH2^2,SH3^3).Sos(dom^3) | (29, 30)

    # Cytosolic
    Shc(PTB,Y317=pY) + Grb2(SH2,SH3) == Shc(PTB,Y317=pY^1).Grb2(SH2^1,SH3) | (31, 32)
    Shc(PTB,Y317=pY) + Grb2(SH2,SH3^_) == Shc(PTB,Y317=pY^1).Grb2(SH2^1,SH3^_) | (33, 34)
    Shc(PTB,Y317=pY) > Shc(PTB,Y317=Y) | 35
    Grb2(SH2,SH3) + Sos(dom) == Grb2(SH2,SH3^1).Sos(dom^1) | (36, 37)
    Shc(PTB,Y317=pY^2).Grb2(SH2^2,SH3) + Sos(dom) == Shc(PTB,Y317=pY^2).Grb2(SH2^2,SH3^3).Sos(dom^3) | (38, 39)


if __name__ == "__main__":
    newseeds = []
    for i, (sp, attr) in enumerate(attributegen()):
        print i, sp, attr
        newseeds.append(sp)
    print ''

    rules = rulegen()
    for i, rr in enumerate(rules):
        print i, rr
    print ''

    generate_reactions(newseeds, rules)
