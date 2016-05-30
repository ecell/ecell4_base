import pythoncyc

class EcocycDataSource(object):
    
    def __init__(self):
        self.pgdb = pythoncyc.select_organism('eco')
        
if __name__ == '__main__':
    eco = EcocycDataSource()
    print(eco.pgdb.all_pathways())
