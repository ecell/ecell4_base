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

def description(entity, collections=None):
    if isinstance(collections, str):
        collections = [collections]

    desc = []

    if collections is None or 'uniprot' in collections:
        from . import uniprot
        desc.extend(uniprot.description(entity))

    if collections is None or 'pdb' in collections:
        from . import pdb
        desc.extend(pdb.description(entity))

    print_descriptions(desc)

def whereis(entity, collections=None):
    if isinstance(collections, str):
        collections = [collections]

    desc = []

    if collections is None or 'uniprot' in collections:
        from . import uniprot
        desc.extend(uniprot.whereis(entity))

    print_descriptions(desc)


__all__ = ['description', 'where']
