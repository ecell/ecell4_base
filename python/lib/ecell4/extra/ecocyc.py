import pythoncyc

class EcocycDataSource(object):
    
    def __init__(self):
        self.pgdb = pythoncyc.select_organism('eco')
        
if __name__ == '__main__':
    from ecell4 import *
    pgdb = EcocycDataSource().pgdb

    with reaction_rules():
        for p in pgdb.all_pathways():
            for g in pgdb.genes_of_pathway(p):
                for p in pgdb.all_products_of_gene(g):
                    ~_eval(g) > _eval(p) | 3

    m = get_model()
    for s in m.list_species():
        print(s.serial())

