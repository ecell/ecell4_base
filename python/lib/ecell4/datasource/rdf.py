from rdflib import Graph, Namespace
from rdflib.namespace import RDF, RDFS, SKOS, DCTERMS, DC
from rdflib.term import URIRef


class RDFDataSourceBase(object):

    GRAPH = {}

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

    def query(self, sparql):
        qres = self.graph.query(sparql)
        for row in qres:
            yield row

    def objects(self, obj, pred):
        qres = self.query(
            """select ?o where
            {{
            ?s
            rdf:type <{:s}>;
            <{:s}> ?o.
            }}
            """.format(obj, pred))
        for row in qres:
            yield row[0]
