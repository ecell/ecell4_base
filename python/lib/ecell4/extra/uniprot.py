from SPARQLWrapper import SPARQLWrapper, JSON

class UniprotDataSource(object):
    def __init__(self):
        self.sparql = SPARQLWrapper("http://sparql.uniprot.org/sparql/")

    def create_species(self, taxon, model):
        self.sparql.setQuery("""

        PREFIX up:<http://purl.uniprot.org/core/> 
        PREFIX taxon:<http://purl.uniprot.org/taxonomy/> 
        PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#> 
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
        SELECT ?protein
        WHERE
        {
        	?protein a up:Protein .
        	?protein up:organism ?organism . 
        	{
        		?protein up:organism taxon:""" + taxon + """.
        	}
        }

        """)
        self.sparql.setReturnFormat(JSON)
        results = self.sparql.query().convert()
        
        for s in results['results']['bindings']:
            uniprot_uri = s['protein']['value']
            uniprot_id = uniprot_uri.split("/")[-1]
            #print(uniprot_id)
            s = Species(uniprot_id)
            s.set_attribute("identifiers.org", "uniprot")
            model.add_species_attribute(s)

if __name__ == '__main__':
    from ecell4 import *
    m1 = NetworkModel()
    ds = UniprotDataSource()
    ds.create_species('83333', m1)

