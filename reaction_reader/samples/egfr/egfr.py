from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.network import generate_reactions


@species_attributes
def attributegen():
    egf_tot = 0.0 
    Grb2_tot= 0.0
    Shc_tot = 0.0
    Sos_tot = 0.0
    egfr_tot= 0.0
    Grb2_Sos_tot = 0.0

    egf(r) | egf_tot
    Grb2(SH2, SH3) | Grb2_tot
    Shc(PTB, Y317=Y) | Shc_tot
    Sos(dom) | Sos_tot
    egfr(l, r, Y1068=Y, Y1148=Y) | egfr_tot
    Grb2(SH2, SH3^1).Sos(dom^1) | Grb2_Sos_tot

@reaction_rules
def rulegen():
    (kp1, km1, kp2, km2, kp3, km3, kp14, km14, kp9, km9, kp11, km11, kp10, km10, kp13, km13, kp15, km15, 
            kp18, km18, kp20, km20, kp17, km17, kp24, km24, kp19, km19,  kp21, km21, kp23, km23, km16, kp12, km12, kp22, km22) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    # Ligand-receptor binding (ligand-monomer)
    egfr(l, r) + egf(r) == egfr(l^1, r).egf(r^1) | (kp1, km1)

    # Note changed multiplicity
    # Receptor-aggregation
    egfr(l^_, r) + egfr(l^_, r) == egfr(l^_,r^3).egfr(l^_,r^3) | (kp2, km2)
    # egfr(l^_1, r) + egfr(l^_2, r) == egfr(l^_1,r^3).egfr(l^_2,r^3) | (kp2, km2) #XXX: this should work, but not now

    # Transphosphorylation of egfr by RTK
    egfr(r^_, Y1068=Y) > egfr(r^_, Y1068=pY) | kp3
    egfr(r^_, Y1148=Y) > egfr(r^_, Y1148=pY) | kp3

    # Dephosphorylation
    egfr(Y1068=pY) > egfr(Y1068=Y) | km3
    egfr(Y1148=pY) > egfr(Y1148=Y) | km3

    # Shc transphosphorylation
    egfr(r^_, Y1148=pY^1).Shc(PTB^1,Y317=Y) > egfr(r^_,Y1148=pY^1).Shc(PTB^1,Y317=pY) | kp14
    Shc(PTB^_,Y317=pY) > Shc(PTB^_,Y317=Y) | km14

    # Y1068 activity
    egfr(Y1068=pY) + Grb2(SH2,SH3) == egfr(Y1068=pY^1).Grb2(SH2^1,SH3) | (kp9, km9)
    egfr(Y1068=pY) + Grb2(SH2,SH3^_) == egfr(Y1068=pY^1).Grb2(SH2^1,SH3^_) | (kp11, km11)
    egfr(Y1068=pY^1).Grb2(SH2^1,SH3) + Sos(dom) == egfr(Y1068=pY^1).Grb2(SH2^1,SH3^2).Sos(dom^2) | (kp10, km10)

    # Y1148 activity
    egfr(Y1148=pY) + Shc(PTB,Y317=Y) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=Y) | (kp13, km13)
    egfr(Y1148=pY) + Shc(PTB,Y317=pY) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY) | (kp15, km15)
    egfr(Y1148=pY) + Shc(PTB,Y317=pY^1).Grb2(SH2^1,SH3) == egfr(Y1148=pY^2).Shc(PTB^2,Y317=pY^1).Grb2(SH2^1,SH3) | (kp18, km18)
    egfr(Y1148=pY) + Shc(PTB,Y317=pY^1).Grb2(SH2^1,SH3^3).Sos(dom^3) == egfr(Y1148=pY^2).Shc(PTB^2,Y317=pY^1).Grb2(SH2^1,SH3^3).Sos(dom^3) | (kp20, km20)

    egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY) + Grb2(SH2,SH3) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY^2).Grb2(SH2^2,SH3) | (kp17, km17)

    egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY) + Grb2(SH2,SH3^3).Sos(dom^3) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY^2).Grb2(SH2^2,SH3^3).Sos(dom^3) | (kp24, km24)

    Shc(PTB^_,Y317=pY^2).Grb2(SH2^2,SH3) + Sos(dom) == Shc(PTB^_,Y317=pY^2).Grb2(SH2^2,SH3^3).Sos(dom^3) | (kp19, km19)

    # Cytosolic
    Shc(PTB,Y317=pY) + Grb2(SH2,SH3) == Shc(PTB,Y317=pY^1).Grb2(SH2^1,SH3) | (kp21, km21)
    Shc(PTB,Y317=pY) + Grb2(SH2,SH3^_) == Shc(PTB,Y317=pY^1).Grb2(SH2^1,SH3^_) | (kp23, km23)
    Shc(PTB,Y317=pY) > Shc(PTB,Y317=Y) | km16
    Grb2(SH2,SH3) + Sos(dom) == Grb2(SH2,SH3^1).Sos(dom^1) | (kp12, km12)
    Shc(PTB,Y317=pY^2).Grb2(SH2^2,SH3) + Sos(dom) == Shc(PTB,Y317=pY^2).Grb2(SH2^2,SH3^3).Sos(dom^3) | (kp22, km22)


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
