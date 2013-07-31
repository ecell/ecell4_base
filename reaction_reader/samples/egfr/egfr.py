import itertools
import ecell4.reaction_reader.species as species

from ecell4.reaction_reader.decorator2 import species_attributes, reaction_rules
from ecell4.reaction_reader.species import FirstOrderReactionRule, SecondOrderReactionRule


def generate_recurse(seeds1, rules, seeds2=[]):
    seeds = list(itertools.chain(seeds1, seeds2))
    retval = []
    for sp1 in seeds1:
        for rr in rules:
            if isinstance(rr, FirstOrderReactionRule):
                try:
                    pttrns = rr.match(sp1)
                except Exception, e:
                    print rr, sp1
                    raise e
                if pttrns is not None and len(pttrns) > 0:
                    for newsp in itertools.chain(*pttrns):
                        if newsp not in seeds and newsp not in retval:
                            retval.append(newsp)
        for sp2 in seeds:
            for rr in rules:
                if isinstance(rr, SecondOrderReactionRule):
                    try:
                        pttrns = rr.match(sp1, sp2)
                    except Exception, e:
                        print rr, sp1, sp2
                        raise e
                    if pttrns is not None and len(pttrns) > 0:
                        for newsp in itertools.chain(*pttrns):
                            if newsp not in seeds and newsp not in retval:
                                retval.append(newsp)
        for sp2 in seeds2:
            for rr in rules:
                if isinstance(rr, SecondOrderReactionRule):
                    try:
                        pttrns = rr.match(sp2, sp1)
                    except Exception, e:
                        print rr, sp1, sp2
                        raise e
                    if pttrns is not None and len(pttrns) > 0:
                        for newsp in itertools.chain(*pttrns):
                            if newsp not in seeds and newsp not in retval:
                                retval.append(newsp)
    return (retval, seeds)

@species_attributes
def attributegen():
    egf(r) | egf_tot
    Grb2(SH2, SH3) | Grb2_tot
    Shc(PTB, Y317=Y) | Shc_tot
    Sos(dom) | Sos_tot
    egfr(l, r, Y1068=Y, Y1148=Y) | egfr_tot
    Grb2(SH2, SH3^1).Sos(dom^1) | Grb2_Sos_tot

@reaction_rules
def rulegen():
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

    # egfr(Y1148=pY^1).Shc(PTB^1,Y317^pY) + Grb2(SH2,SH3^3).Sos(dom^3) == egfr(Y1148=pY^1).Shc(PTB^1,Y317=pY^2).Grb2(SH2^2,SH3^3).Sos(dom^3) | (kp24, km24)
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

    seeds, cnt = [], 0
    while len(newseeds) != 0 and cnt < 10:
        print "[RESULT%d: %d]" % (cnt, len(seeds)), newseeds, seeds
        newseeds, seeds = generate_recurse(newseeds, rules, seeds)
        cnt += 1
    print "[RESULT%d: %d]" % (cnt, len(seeds)), newseeds, seeds
