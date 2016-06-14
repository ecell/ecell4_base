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
PREFIX skos:<http://www.w3.org/2004/02/skos/core#>
SELECT ?protein ?genename
WHERE
{
    ?protein a up:Protein .
    ?protein up:organism ?organism .
    ?protein up:organism taxon:""" + taxon + """ .
    ?protein up:encodedBy ?gene .
    ?gene skos:prefLabel ?genename .
}

        """)
        self.sparql.setReturnFormat(JSON)
        results = self.sparql.query().convert()
        
        for s in results['results']['bindings']:
            uniprot_uri = s['protein']['value']
            uniprot_id = uniprot_uri.split("/")[-1]
            uniprot_genename = s['genename']['value']
            s = Species(uniprot_id)
#            if s.serial() == uniprot_id:
#                print(uniprot_id)
#            else:
            s.set_attribute("identifiers.org", "uniprot")
            s.set_attribute("genename", uniprot_genename)
            model.add_species_attribute(s)

if __name__ == '__main__':
    from ecell4 import *
    m1 = NetworkModel()
    ds = UniprotDataSource()
    ds.create_species('4932', m1)

