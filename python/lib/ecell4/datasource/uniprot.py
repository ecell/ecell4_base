import re

from rdflib import Namespace
from rdflib.namespace import RDF, RDFS, SKOS
from rdflib.term import URIRef

try:
    from . import rdf
    from . import pdb
except SystemError:
    import rdf
    import pdb

class UniProtDataSourceBase(rdf.RDFDataSourceBase):

    GRAPH = {}
    UNIPROT = Namespace("http://purl.uniprot.org/core/")
    UPDB = Namespace("http://purl.uniprot.org/database/")

    def __init__(self, url=None, cache=True):
        rdf.RDFDataSourceBase.__init__(self, url, cache)

class UniProtTaxonDataSource(UniProtDataSourceBase):

    URL = "http://www.uniprot.org/taxonomy/{entry_id}.rdf"

    def __init__(self, entry_id=None, url=None, cache=True):
        if url is not None:
            UniProtDataSourceBase.__init__(
                self, url, cache)
        elif entry_id is not None:
            UniProtDataSourceBase.__init__(
                self, self.URL.format(entry_id=entry_id), cache)
        else:
            UniProtDataSourceBase.__init__(self, None, cache)

    def scientific_name(self):
        return [str(obj) for obj in self.graph.objects(predicate=self.UNIPROT.scientificName)]

class UniProtDataSource(UniProtDataSourceBase):

    URL = "http://www.uniprot.org/uniprot/{entry_id}.rdf"

    def __init__(self, entry_id=None, url=None, cache=True):
        if url is not None:
            UniProtDataSourceBase.__init__(
                self, url, cache)
        elif entry_id is not None:
            UniProtDataSourceBase.__init__(
                self, self.URL.format(entry_id=entry_id), cache)
        else:
            UniProtDataSourceBase.__init__(self, None, cache)

    def gene(self):
        return [str(obj) for obj in self.objects(self.UNIPROT.Gene, SKOS.prefLabel)]

    def locus_name(self):
        return [str(obj) for obj in self.objects(self.UNIPROT.Gene, self.UNIPROT.locusName)]

    def function_annotation(self):
        return [str(obj) for obj in self.objects(self.UNIPROT.Function_Annotation, RDFS.comment)]

    def organism(self):
        retval = []
        for obj1 in self.graph.objects(predicate=self.UNIPROT.organism):
            mobj = re.match("http:\/\/purl\.uniprot\.org\/taxonomy\/([0-9]+)", str(obj1))
            if mobj is None:
                continue
            # taxon_id = mobj.group(1)
            # retval.extend(UniProtTaxonDataSource(taxon_id).scientific_name())
            retval.extend(UniProtTaxonDataSource(url=str(obj1)).scientific_name())
        return retval

    def structure_resource(self):
        return [str(sub) for sub in self.graph.subjects(predicate=RDF.type, object=self.UNIPROT.Structure_Resource)]

    def pdb(self):
        retval = []
        for sub in self.graph.subjects(predicate=RDF.type, object=self.UNIPROT.Structure_Resource):
            if URIRef("http://purl.uniprot.org/database/PDB") not in self.graph.objects(subject=sub, predicate=self.UNIPROT.database):
                continue
            # mobj = re.match("http:\/\/rdf\.wwpdb\.org\/pdb\/([0-9A-Za-z]+)", str(sub))
            # if mobj is None:
            #     continue
            # pdb_id = mobj.group(1).upper()
            retval.extend(pdb.PDBDataSource(url=str(sub)).identifier())
        return retval

    def database(self, dbname):
        return [str(sub) for sub in self.graph.subjects(predicate=self.UNIPROT.database, object=self.UPDB[dbname])]

    def biogrid(self):
        return self.database("BioGrid")


if __name__ == "__main__":
    print(UniProtDataSource("P0AEZ3").gene()[0])
    print(UniProtDataSource("P0AEZ3").locus_name())
    print(UniProtDataSource("P0AEZ3").function_annotation()[0])
    print(UniProtDataSource("P0AEZ3").organism()[0])
    print(UniProtDataSource("P0AEZ3").structure_resource())
    print(UniProtDataSource("P0AEZ3").pdb())
    print(UniProtDataSource("P0AEZ3").biogrid())
    print(UniProtDataSource("P0AEZ3").database("IntAct"))
