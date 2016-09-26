def description(entity, collections=None):
    if isinstance(collections, str):
        collections = [collections]

    if collections is None or 'uniprot' in collections:
        from . import uniprot
        for line in uniprot.description(entity):
            print(line)

    if collections is None or 'pdb' in collections:
        from . import pdb
        for line in pdb.description(entity):
            print(line)


__all__ = ['description']
