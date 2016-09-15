import re

from rdflib import Graph, Namespace
from rdflib.namespace import RDF, RDFS, SKOS, DCTERMS, DC
from rdflib.term import URIRef


class PDBDataSourceBase(object):

    GRAPH = {}
    UNIPROT = Namespace("http://purl.uniprot.org/core/")

    def __init__(self, url=None, cache=True):
        self.url = url
        self.cache = cache

        if self.url is not None:
            self.graph = self.fetch(self.url, self.cache)

    def fetch(self, url, cache=False):
        if not cache or url not in self.GRAPH.keys():
            graph = Graph()
            graph.parse(url, format="xml")

            if cache:
                self.GRAPH[url] = graph
        else:
            graph = self.GRAPH[url]
        return graph

    def objects(self, uri, pred):
        for sub in self.graph.subjects(predicate=RDF.type, object=uri):
            for obj in self.graph.objects(subject=sub, predicate=pred):
                yield obj

class PDBDataSource(PDBDataSourceBase):

    URL = "http://rdf.wwpdb.org/pdb/{entry_id}"

    def __init__(self, entry_id=None, url=None, cache=True):
        if url is not None:
            PDBDataSourceBase.__init__(
                self, url, cache)
        elif entry_id is not None:
            PDBDataSourceBase.__init__(
                self, self.URL.format(entry_id=entry_id), cache)
        else:
            PDBDataSourceBase.__init__(self, None, cache)

    def identifier(self):
        return [str(obj) for obj in self.graph.objects(predicate=DCTERMS.identifier)]

    def title(self):
        return [str(obj) for obj in self.graph.objects(predicate=DC.title)]


if __name__ == "__main__":
    print(PDBDataSource("3Q9L").identifier())
    print(PDBDataSource("3Q9L").title())
