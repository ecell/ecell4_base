import re

from rdflib import Namespace
from rdflib.namespace import DCTERMS, DC

try:
    from . import rdf
except SystemError:
    import rdf


class PDBDataSource(rdf.RDFDataSourceBase):

    UNIPROT = Namespace("http://purl.uniprot.org/core/")
    URL = "http://rdf.wwpdb.org/pdb/{entry_id}"

    def __init__(self, entry_id=None, url=None, cache=True):
        if url is not None:
            rdf.RDFDataSourceBase.__init__(
                self, url, cache)
        elif entry_id is not None:
            rdf.RDFDataSourceBase.__init__(
                self, self.URL.format(entry_id=entry_id), cache)
        else:
            rdf.RDFDataSourceBase.__init__(self, None, cache)

    def identifier(self):
        return [str(obj) for obj in self.graph.objects(predicate=DCTERMS.identifier)]

    def title(self):
        return [str(obj) for obj in self.graph.objects(predicate=DC.title)]


if __name__ == "__main__":
    print(PDBDataSource("3Q9L").identifier())
    print(PDBDataSource("3Q9L").title())
