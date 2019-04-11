from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
#from ecell4.reaction_reader.species import generate_reactions
from ecell4.util.sbml_exporter import convert2SBML
from ecell4.reaction_reader.network import generate_NetworkModel, generate_reactions


@species_attributes
def attributegen():
    egf(r) | 0.1
    Grb2(SH2, SH3) | 0.2
    Shc(PTB, Y317=Y) | 0.3
    Sos(dom) | 0.4
    egfr(l, r, Y1068=Y, Y1148=Y) | 0.5
    Grb2(SH2, SH3^1).Sos(dom^1) | 0.6

@reaction_rules
def rulegen():
    # Ligand-receptor binding (ligand-monomer)
    egfr(l, r) + egf(r) == egfr(l^1, r).egf(r^1) | (0.1, 0.2)

    # Note changed multiplicity
    # Receptor-aggregation
    egfr(l^_, r) + egfr(l^_, r) == egfr(l^_,r^3).egfr(l^_,r^3) | (0.3, 0.4)
    # egfr(l^_1, r) + egfr(l^_2, r) == egfr(l^_1,r^3).egfr(l^_2,r^3) | (kp2, km2) #XXX: this should work, but not now

    # Transphosphorylation of egfr by RTK
    egfr(r^_, Y1068=Y) > egfr(r^_, Y1068=pY) | 0.5
    egfr(r^_, Y1148=Y) > egfr(r^_, Y1148=pY) | 0.6

    # Dephosphorylation
    egfr(Y1068=pY) > egfr(Y1068=Y) | 0.7
    egfr(Y1148=pY) > egfr(Y1148=Y) | 0.8

    # Shc transphosphorylation
    egfr(r^_, Y1148=pY^1).Shc(PTB^1,Y317=Y) > egfr(r^_,Y1148=pY^1).Shc(PTB^1,Y317=pY) | 0.9
    Shc(PTB^_,Y317=pY) > Shc(PTB^_,Y317=Y) | 1.0

    # Y1068 activity
    egfr(Y1068=pY) + Grb2(SH2,SH3) == egfr(Y1068=pY^1).Grb2(SH2^1,SH3) | (0.1, 0.2)
    egfr(Y1068=pY) + Grb2(SH2,SH3^_) == egfr(Y1068=pY^1).Grb2(SH2^1,SH3^_) | (0.3, 0.4)
    egfr(Y1068=pY^1).Grb2(SH2^1,SH3) + Sos(dom) == egfr(Y1068=pY^1).Grb2(SH2^1,SH3^2).Sos(dom^2) | (0.5, 0.6)

    # Y1148 activity
    egfr(Y1148=pY) + Shc(PTB,Y317=Y) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=Y) | (0.7, 0.8)
    egfr(Y1148=pY) + Shc(PTB,Y317=pY) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY) | (0.9, 1.0)
    egfr(Y1148=pY) + Shc(PTB,Y317=pY^1).Grb2(SH2^1,SH3) == egfr(Y1148=pY^2).Shc(PTB^2,Y317=pY^1).Grb2(SH2^1,SH3) | (0.1, 0.2)
    egfr(Y1148=pY) + Shc(PTB,Y317=pY^1).Grb2(SH2^1,SH3^3).Sos(dom^3) == egfr(Y1148=pY^2).Shc(PTB^2,Y317=pY^1).Grb2(SH2^1,SH3^3).Sos(dom^3) | (0.3, 0.4)

    egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY) + Grb2(SH2,SH3) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY^2).Grb2(SH2^2,SH3) | (0.5, 0.6)

    egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY) + Grb2(SH2,SH3^3).Sos(dom^3) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY^2).Grb2(SH2^2,SH3^3).Sos(dom^3) | (0.7, 0.8)

    Shc(PTB^_,Y317=pY^2).Grb2(SH2^2,SH3) + Sos(dom) == Shc(PTB^_,Y317=pY^2).Grb2(SH2^2,SH3^3).Sos(dom^3) | (0.9, 1.0)

    # Cytosolic
    Shc(PTB,Y317=pY) + Grb2(SH2,SH3) == Shc(PTB,Y317=pY^1).Grb2(SH2^1,SH3) | (1.1, 1.2)
    Shc(PTB,Y317=pY) + Grb2(SH2,SH3^_) == Shc(PTB,Y317=pY^1).Grb2(SH2^1,SH3^_) | (1.3, 1.4)
    Shc(PTB,Y317=pY) > Shc(PTB,Y317=Y) | 1.5
    Grb2(SH2,SH3) + Sos(dom) == Grb2(SH2,SH3^1).Sos(dom^1) | (1.6, 1.7)
    Shc(PTB,Y317=pY^2).Grb2(SH2^2,SH3) + Sos(dom) == Shc(PTB,Y317=pY^2).Grb2(SH2^2,SH3^3).Sos(dom^3) | (1.8, 1.9)


if __name__ == "__main__":
    newseeds = []
    attrs = attributegen()
    for i, (sp, attr) in enumerate(attrs):
        #print i, sp, attr
        newseeds.append(sp)
        #print ''
    reaction_rules = rulegen()

    seeds, rules = generate_reactions(newseeds , reaction_rules, max_iter = 5)

    m = generate_NetworkModel(seeds, rules)
    convert2SBML(m, attrs, "egfr_n5.xml")
