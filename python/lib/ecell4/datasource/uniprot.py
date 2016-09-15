import re

# try:
#     from urllib2 import Request, urlopen
# except ImportError:
#     from urllib.request import Request, urlopen

from rdflib import Graph, Namespace
from rdflib.namespace import RDF, RDFS, SKOS
from rdflib.term import URIRef

try:
    from . import pdb
except SystemError:
    import pdb

class UniProtDataSourceBase(object):

    GRAPH = {}
    UNIPROT = Namespace("http://purl.uniprot.org/core/")
    UPDB = Namespace("http://purl.uniprot.org/database/")

    def __init__(self, url=None, cache=True):
        self.url = url
        self.cache = cache

        if self.url is not None:
            self.graph = self.fetch(self.url, self.cache)

    def fetch(self, url, cache=False):
        if not cache or url not in self.GRAPH.keys():
            # req = Request(url)
            # response = urlopen(req)
            # rdf = response.read().decode('utf-8')
            graph = Graph()
            graph.parse(url)

            if cache:
                self.GRAPH[url] = graph
        else:
            graph = self.GRAPH[url]
        # for s, p, o in graph:
        #     print(repr(s), repr(p), repr(o))
        return graph

    def objects(self, uri, pred):
        # for sub in self.graph.subjects(predicate=RDF.type, object=uri):
        #     for obj in self.graph.objects(subject=sub, predicate=pred):
        #         yield obj
        qres = self.graph.query(
            """select ?obj where
            {{
            ?s
            rdf:type <{:s}>;
            <{:s}> ?obj.
            }}
            """.format(uri, pred))
        for row in qres:
            yield row[0]

class UniProtTaxonDataSource(UniProtDataSourceBase):

    URL = "http://www.uniprot.org/taxonomy/{entry_id}.rdf"

    def __init__(self, entry_id=None, cache=True):
        if entry_id is not None:
            UniProtDataSourceBase.__init__(
                self, self.URL.format(entry_id=entry_id), cache)
            self.entry_id = entry_id
        else:
            UniProtDataSourceBase.__init__(self, None, cache)
            self.entry_id = None

    def scientific_name(self):
        return [str(obj) for obj in self.graph.objects(predicate=self.UNIPROT.scientificName)]

class UniProtDataSource(UniProtDataSourceBase):

    URL = "http://www.uniprot.org/uniprot/{entry_id}.rdf"

    def __init__(self, entry_id=None, cache=True):
        if entry_id is not None:
            UniProtDataSourceBase.__init__(
                self, self.URL.format(entry_id=entry_id), cache)
            self.entry_id = entry_id
        else:
            UniProtDataSourceBase.__init__(self, None, cache)
            self.entry_id = None

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
            taxon_id = mobj.group(1)
            retval.extend(UniProtTaxonDataSource(taxon_id).scientific_name())
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
