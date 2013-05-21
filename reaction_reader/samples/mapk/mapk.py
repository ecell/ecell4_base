from ecell4.reaction_reader.decorators import reaction_rules


@reaction_rules
def reactions(kon, koff, kcat):
    mapk(phos=YT) + kk(bs=on) > mapk(phos=YT[1]).kk(bs=on[1]) | kon
    mapk(phos=YT[1]).kk(bs=on[1]) > mapk(phos=YT) + kk(bs=on) | koff
    mapk(phos=YT[1]).kk(bs=on[1]) > mapk(phos=pYT) + kk(bs=on) | kcat

    mapk(phos=pYT) + pp(bs=on) > mapk(phos=pYT[1]).pp(bs=on[1]) | kon
    mapk(phos=pYT[1]).pp(bs=on[1]) > mapk(phos=pYT) + pp(bs=on) | koff
    mapk(phos=pYT[1]).pp(bs=on[1]) > mapk(phos=YT) + pp(bs=on) | kcat


if __name__ == "__main__":
    rules = reactions(1, 2, 3)
    print rules
