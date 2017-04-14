def print_descriptions(desc):
    for i, entry in enumerate(desc):
        if i > 0:
            print()
        for line in entry:
            assert len(line) == 2 or len(line) == 3
            if line[1] is None or line[1] == '':
                continue
            if len(line) == 2:
                print('{0}: {1}'.format(*line))
            else:
                print('{0}{2}{1}'.format(*line))

def __description(entity, collections):
    desc = []

    if collections is None or 'uniprot' in collections:
        from . import uniprot
        desc.extend(uniprot.description(entity))

    if collections is None or 'pdb' in collections:
        from . import pdb
        desc.extend(pdb.description(entity))

    if collections is None or 'pubmed' in collections:
        from . import pubmed
        desc.extend(pubmed.description(entity))

    return desc

def description(entity, collections=None):
    from ecell4 import Species

    if isinstance(collections, str):
        collections = [collections]

    if isinstance(entity, (str, Species)):
        desc = __description(entity, collections)
    else:
        desc = []
        for e in entity:
            desc.extend(__description(e, collections))

    print_descriptions(desc)

def whereis(entity, collections=None):
    if isinstance(collections, str):
        collections = [collections]

    desc = []

    if collections is None or 'uniprot' in collections:
        from . import uniprot
        desc.extend(uniprot.whereis(entity))

    print_descriptions(desc)


__all__ = ['description', 'whereis']
