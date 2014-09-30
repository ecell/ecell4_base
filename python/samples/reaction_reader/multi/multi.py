from ecell4.reaction_reader.decorator2 import reaction_rules, create_species
from ecell4.reaction_reader.network import generate_reactions


def test1():
    @reaction_rules
    def rulegen():
        _(ps=u) > _(ps=p)
        _1(bs) + _1(bs) > _1(bs^1)._1(bs^1)
        # (_(ps=u) + K(bs) == _(ps=u^1).K(bs^1) | ExcludeReactants(1, K)
        #     > _(ps=p) + K(bs))

    rules = rulegen()
    print generate_reactions(
        [create_species("K(bs)"), create_species("A(ps=u)")], rules)[0]
    print generate_reactions(
        [create_species("K(bs)"), create_species("A(ps1=u,ps2=u)")], rules)[0]
    print generate_reactions(
        [create_species("K(bs)"),
        create_species("A(ps1=u,ps2=u,ps=(ps1,ps2))")], rules)[0]
    print generate_reactions(
        [create_species("K(bs)"),
        create_species("A(ps1=u,ps2=u,ps=(ps1,ps2),bs=(ps1,))")], rules)[0]

def test2():
    @reaction_rules
    def rulegen():
        (A(bs) + B(bs) > A(bs^1).B(bs^1)
            | ExcludeReactants(1, B) | ExcludeReactants(2, A))

    rules = rulegen()
    print generate_reactions(
        [create_species("A(bs1, bs2, bs=(bs1, bs2))"),
        create_species("B(bs1, bs3, bs=(bs1, bs3))")], rules)[0]

def test3():
    @reaction_rules
    def rulegen():
        # _(ps=u) > _(ps=p)
        # _(_1=u) > _(_1=p)
        # A(ps1=u) > A(ps1=p)
        # A(ps1=u,ps=(ps1,ps2)) > A(ps1=p,ps=(ps1,ps2))
        # A(_1=u,ps=(_1,_)) > A(_1=p,ps=(_1,_))
        A(_1=u,_2=u,ps=(_1,_2)) > A(_1=p,_2=u,ps=(_1,_2))

    rules = rulegen()
    print generate_reactions(
        [create_species("A(ps1=u,ps=(ps1,))")], rules)[0]
    print generate_reactions(
        [create_species("A(ps1=u,ps2=u,ps=(ps1,ps2))")], rules)[0]
    print generate_reactions(
        [create_species("A(ps1=u,ps2=u,ps3=u,ps=(ps1,ps2,ps3))")], rules)[0]

def test4():
    @reaction_rules
    def rulegen():
        A(_1=u,ps=(_1,)) + A(_1=u,ps=(_1,)) > A(_1=u^1,ps=(_1,)).A(_1=u^1,ps=(_1,))

    rules = rulegen()
    print generate_reactions(
        [create_species("A(ps1=u,ps2=u,ps3=u,ps=(ps1,ps2,ps3))")], rules, max_iter=1)[0]


if __name__ == '__main__':
    # test1()
    # test2()
    # test3()
    test4()
