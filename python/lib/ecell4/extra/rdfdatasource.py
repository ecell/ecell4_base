from SPARQLWrapper import SPARQLWrapper, JSON
from collections import defaultdict

class ReactomeDataSource(object):
    def __init__(self):
        self.sparql = SPARQLWrapper("https://www.ebi.ac.uk/rdf/services/reactome/sparql")

    def create_reactions(self, taxon, model):
        self.sparql.setQuery("""

PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX dc: <http://purl.org/dc/elements/1.1/>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX foaf: <http://xmlns.com/foaf/0.1/>
PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX biopax3: <http://www.biopax.org/release/biopax-level3.owl#>

SELECT DISTINCT ?pathway ?pathwayname ?rea ?lid ?rid
WHERE
{
 ?pathway rdf:type biopax3:Pathway . 
 ?pathway biopax3:displayName ?pathwayname . 
 ?pathway biopax3:organism <http://identifiers.org/taxonomy/""" + taxon + """> .
 ?pathway biopax3:pathwayComponent ?rea .
 ?rea biopax3:left ?l .
 ?l biopax3:entityReference ?lid .
 ?rea biopax3:right ?r .
 ?r biopax3:entityReference ?rid .
}

        """)
        self.sparql.setReturnFormat(JSON)
        results = self.sparql.query().convert()

        lefts = defaultdict(list)
        rights = defaultdict(list)
        for r in results['results']['bindings']:
            lefts[r['rea']['value']].append(r['lid']['value'])
            rights[r['rea']['value']].append(r['rid']['value'])

        rr = ReactionRule()
        for r in lefts.keys():
            for l in set(lefts[r]):
                #print(l.split("/")[-1])
                try:
                    rr.add_reactant(Species(l.split("/")[-1]))
                except:
                    pass
            for r in set(rights[r]):
                #print(r.split("/")[-1])
                try:
                    rr.add_product(Species(r.split("/")[-1]))
                except:
                    pass
        model.add_reaction_rule(rr)

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
    upds = UniprotDataSource()
    upds.create_species('4932', m1)
    rtds = ReactomeDataSource()
    rtds.create_reactions('4932', m1)


