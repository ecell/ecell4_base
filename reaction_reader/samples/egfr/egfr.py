from ecell4.reaction_reader.decorator import just_parse, reaction_rules

@just_parse
def reactions(kp1, km1, kp2, km2, kp3, km3):
    # Ligand-receptor binding (ligand-monomer)
    egfr(l,r) + egf(r) == egfr(l^1,r).egf(r^1) | (kp1, km1)

    # Note changed multiplicity
    # Receptor-aggregation
    egfr(l^_1, r) + egfr(l^_2, r) == egfr(l^_1,r^3).egfr(l^_2,r^3) | (kp2, km2)

    # Transphosphorylation of egfr by RTK
    egfr(r^_, Y1068=Y) > egfr(r^_, Y1068=pY) | kp3
    egfr(r^_, Y1148=Y) > egfr(r^_, Y1148=pY) | kp3

    # Dephosphorylation
    egfr(Y1068=pY) > egfr(Y1068=Y) | km3
    egfr(Y1148=pY) > egfr(Y1148=Y) | km3

    # Shc transphosphorylation
    egfr(r^_, Y1148=pY^1).Shc(PTB^1,Y317=Y) > egfr(r^_,Y1148=pY^1).Shc(PTB^1,Y317=pY) | kp14
    Shc(PTB^_,Y317=pY) > Shc(PTB^_,Y317=Y) | km14

if __name__ == "__main__":
    rules = reactions(1, 2, 3, 4, 5, 6)
    for i, rr in enumerate(rules):
        print i + 1, rr
