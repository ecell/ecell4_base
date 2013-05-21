from ecell4.reaction_reader.decorators import reaction_rules


@reaction_rules
def reactions(kf, kr):
    mapk(phos=YT) + kk(bs=on) > mapk(phos=YT[1]).kk(bs=on[1]) | kf
    mapk(phos=YT[1]).kk(bs=on[1]) > mapk(phos=pYT) + kk(bs=on) | kr


if __name__ == "__main__":
    rules = reactions(1, 2)
    print rules
