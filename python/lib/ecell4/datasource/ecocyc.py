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

        # for x in pgdb.reactions.instances:
        #     left = pgdb.get_slot_value(x['frameid'], 'Left')
        #     right = pgdb.get_slot_value(x['frameid'], 'Right')
        #     if left == None or right == None:
        #         print x['frameid']
        #     else:
        #         left = left.strip('|').strip()
        #         right = right.strip('|').strip()
        #         if "-" in left:
        #             left = left.replace("-", "__")
        #         if "-" in right:
        #             right = right.replace("-", "__")
        #         if "+" in left:
        #             left = left.replace("+", "__")
        #         if "+" in right:
        #             right = right.replace("+", "__")
        #         if re.search('^[0-9]', left):
        #             left = '__' + left
        #         if re.search('^[0-9]', right):
        #             right = '__' + right

        #         _eval(left) > _eval(right) | 10

    m = get_model()
    for s in m.list_species():
        print(s.serial())

