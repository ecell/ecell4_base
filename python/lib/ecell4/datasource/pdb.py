import re

from rdflib import Namespace
from rdflib.namespace import DCTERMS, DC, RDFS

try:
    from . import rdf
except SystemError:
    import rdf

from rdflib import Graph, Namespace

def description(entity):
    entity_id = PDBDataSource.parse_entity(entity)
    if entity_id is not None:
        entry = []
        src = PDBDataSource(entity_id)
        entry.append(("PDB", ', '.join(src.identifier()), ' - '))
        entry.append(("Title", ', '.join(src.title())))
        # desc.append("Protein: {}".format(', '.join(src.structured_name())))
        src_gen = src.src_gen()
        if len(src_gen) > 0:
            entry.append(("Gene", src_gen[0]["gene"]))
            entry.append(("Organism", "{} {}".format(src_gen[0]["scientific_name"], src_gen[0]["strain"])))
        for url in src.see_also():
            entry.append(("See Also", url))
        entry.append(("URL", PDBDataSource.link(entity)))
        return [entry]
    return []

class PDBDataSource(rdf.RDFDataSourceBase):

    PDBo = Namespace("http://rdf.wwpdb.org/schema/pdbx-v40.owl#")
    URL = "http://rdf.wwpdb.org/pdb/{entity_id}"

    def __init__(self, entity=None, url=None, cache=True):
        if url is not None:
            rdf.RDFDataSourceBase.__init__(
                self, url, cache)
        elif entity is not None:
            entity_id = self.parse_entity(entity)
            if entity_id is not None:
                rdf.RDFDataSourceBase.__init__(
                    self, self.URL.format(entity_id=entity_id), cache)
            else:
                rdf.RDFDataSourceBase.__init__(self, None, cache)
        else:
            rdf.RDFDataSourceBase.__init__(self, None, cache)

    @classmethod
    def parse_entity(cls, entity):
        # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000020
        collection = 'pdb'
        idpttrn = r'[0-9][A-Za-z0-9]{3}'
        uri1 = r'http://rdf.wwpdb.org/pdb/(?P<id>{})'.format(idpttrn)
        if isinstance(entity, str):
            if re.match(r'^{}$'.format(idpttrn), entity) is not None:
                return entity
            mobj = re.match(uri1, entity)
            if mobj is not None:
                return mobj.group('id')
        else:
            import ecell4
            if isinstance(entity, ecell4.Species) and entity.has_attribute(collection):
                return cls.parse_entity(entity.get_attribute(collection))
        return None  #XXX: Error

    @classmethod
    def link(cls, entity):
        entity_id = cls.parse_entity(entity)
        assert entity_id is not None
        return "http://rdf.wwpdb.org/pdb/{}".format(entity_id)

    def identifier(self):
        return [str(obj) for obj in self.graph.objects(predicate=DCTERMS.identifier)]

    def title(self):
        return [str(obj) for obj in self.graph.objects(predicate=DC.title)]

    def see_also(self):
        return [str(obj) for obj in self.graph.objects(predicate=RDFS.seeAlso)]

    def src_gen(self):
        retval = []
        for obj in self.graph.objects(predicate=self.PDBo.has_entity_src_genCategory):
            qres = rdf.RDFDataSourceBase(url=str(obj), cache=self.cache).query(
                """prefix PDBo: <http://rdf.wwpdb.org/schema/pdbx-v40.owl#>
                select ?entity_id ?scientific_name ?strain ?gene where
                {{
                ?s PDBo:has_entity_src_gen ?entity .
                ?entity
                PDBo:entity_src_gen.entity_id ?entity_id ;
                PDBo:entity_src_gen.pdbx_gene_src_scientific_name ?scientific_name ;
                PDBo:entity_src_gen.gene_src_strain ?strain ;
                PDBo:entity_src_gen.pdbx_gene_src_gene ?gene .
                }}
                """)
            for row in qres:
                entity_id, scientific_name, strain, gene = [str(elem) for elem in row]
                retval.append(
                    dict(entity_id=entity_id, scientific_name=scientific_name, strain=strain, gene=gene))
        return retval


if __name__ == "__main__":
    print(PDBDataSource("3Q9L").identifier())
    print(PDBDataSource("3Q9L").title())
    print(PDBDataSource("3Q9L").see_also())
    print(PDBDataSource("3Q9L").src_gen())

    print(description("3Q9L"))
