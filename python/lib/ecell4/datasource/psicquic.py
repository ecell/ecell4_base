from rdflib.namespace import RDF
from rdflib import Namespace

try:
    from . import rdf
except SystemError:
    import rdf


class PSICQUICDataSource(rdf.RDFDataSourceBase):

    BIOPAX = Namespace("http://www.biopax.org/release/biopax-level3.owl#")
    URL = "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/{entry_id}?format=rdf-xml"

    def __init__(self, entry_id=None, url=None, cache=True):
        if url is not None:
            rdf.RDFDataSourceBase.__init__(
                self, url, cache)
        elif entry_id is not None:
            rdf.RDFDataSourceBase.__init__(
                self, self.URL.format(entry_id=entry_id), cache)
        else:
            rdf.RDFDataSourceBase.__init__(self, None, cache)

    def protein(self):
        return [str(sub) for sub in self.graph.subjects(RDF.type, self.BIOPAX.ProteinReference)]

    def small_molecule(self):
        return [str(sub) for sub in self.graph.subjects(RDF.type, self.BIOPAX.SmallMoleculeReference)]

    def interactor(self):
        return self.protein() + self.small_molecule()

    def interaction(self):
        return [str(sub) for sub in self.graph.subjects(RDF.type, self.BIOPAX.MolecularInteraction)]


if __name__ == "__main__":
    res = PSICQUICDataSource("P0AEZ3").protein()
    print(res, len(res))
    res = PSICQUICDataSource("P0AEZ3").small_molecule()
    print(res, len(res))
    res = PSICQUICDataSource("P0AEZ3").interaction()
    print(res, len(res))
