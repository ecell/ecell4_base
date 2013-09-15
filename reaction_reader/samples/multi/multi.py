from ecell4.reaction_reader.decorator2 import reaction_rules, create_species
from ecell4.reaction_reader.network import generate_reactions


@reaction_rules
def rulegen():
    _(ps=u) > _(ps=p)
    # A(ps=_1) > A(ps=_1)


if __name__ == '__main__':
    rules = rulegen()

    print generate_reactions([create_species("A(ps=u)")], rules)
    print generate_reactions([create_species("A(ps1=u, ps2=u)")], rules)
    print generate_reactions([create_species("A(ps1=u, ps2=u, ps=(ps1, ps2))")], rules)
    species, reactions = generate_reactions(
        [create_species("A(ps1=u, ps2=u, ps3=u, ps=(ps1, ps2, ps3))")], rules)
    print (species, reactions)
    print len(species)
