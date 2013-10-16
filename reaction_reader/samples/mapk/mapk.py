from ecell4.reaction_reader.decorator_base import just_parse
from ecell4.reaction_reader.decorator import reaction_rules


# @reaction_rules
@just_parse
def reactions(kon, koff, kcat):
    mapk(phos=YT) + kk(bs) > mapk(phos=YT^1).kk(bs^1) | kon
    mapk(phos=YT^1).kk(bs^1) > mapk(phos=YT) + kk(bs) | koff
    mapk(phos=YT^1).kk(bs^1) > mapk(phos=pYT) + kk(bs) | kcat

    mapk(phos=pYT) + pp(bs) == mapk(phos=pYT^1).pp(bs^1) | (kon, koff)
    mapk(phos=pYT^1).pp(bs^1) > mapk(phos=YT) + pp(bs) | kcat

    mapk(phos=pYT) + kk(bs) <> mapk(phos=pYT^1).kk(bs^1) | (kon, koff)
    mapk(phos=pYT^1).kk(bs^1) > mapk(phos=pYpT) + kk(bs) | kcat

    (mapk(phos=pYpT) + pp(bs)
        == mapk(phos=pYpT^1).pp(bs^1) | (kon, koff)
        > mapk(phos=pYT) + pp(bs) | kcat)

    # (mapk(phos=YT) + kk(bs)
    #     == mapk(phos=YT^1).kk(bs^1) | (kon, koff)
    #     > mapk(phos=pYT) + kk(bs) | kcat
    #     == mapk(phos=pYT^1).kk(bs^1) | (kon, koff)
    #     > mapk(phos=pYpT) + kk(bs) | kcat)

    # (mapk(phos=pYpT) + pp(bs)
    #     == mapk(phos=pYpT^1).pp(bs^1) | (kon, koff)
    #     > mapk(phos=pYT) + pp(bs) | kcat
    #     == mapk(phos=pYT^1).pp(bs^1) | (kon, koff)
    #     > mapk(phos=YT) + pp(bs) | kcat)


if __name__ == "__main__":
    rules = reactions(1, 2, 3)
    for i, rr in enumerate(rules):
        print i + 1, rr
