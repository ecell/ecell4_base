from ecell4.reaction_reader.decorator import reaction_rules


@reaction_rules
def reactions(kon, koff, kcat):
    mapk(phos=YT) + kk(bs) > mapk(phos=YT[1]).kk(bs[1]) | kon
    mapk(phos=YT[1]).kk(bs[1]) > mapk(phos=YT) + kk(bs) | koff
    mapk(phos=YT[1]).kk(bs[1]) > mapk(phos=pYT) + kk(bs) | kcat

    mapk(phos=pYT) + pp(bs) > mapk(phos=pYT[1]).pp(bs[1]) | kon
    mapk(phos=pYT[1]).pp(bs[1]) > mapk(phos=pYT) + pp(bs) | koff
    mapk(phos=pYT[1]).pp(bs[1]) > mapk(phos=YT) + pp(bs) | kcat


if __name__ == "__main__":
    rules = reactions(1, 2, 3)
    for i, rr in enumerate(rules):
        # print i + 1, rr
        print i + 1, (rr[0] >> list()), (rr[1] >> list())
