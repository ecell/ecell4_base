try:
    from . import rdf
except SystemError:
    import rdf


class PSICQUICDataSource(rdf.RDFDataSourceBase):

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
        qres =  self.query(
             """prefix biopax: <http://www.biopax.org/release/biopax-level3.owl#>
             select ?s where {{ ?s rdf:type biopax:ProteinReference . }}
             """)
        return [str(row[0]) for row in qres]

    def small_molecule(self):
        qres =  self.query(
             """prefix biopax: <http://www.biopax.org/release/biopax-level3.owl#>
             select ?s where {{ ?s rdf:type biopax:SmallMoleculeReference . }}
             """)
        return [str(row[0]) for row in qres]

    def interaction(self):
        qres =  self.query(
             """prefix biopax: <http://www.biopax.org/release/biopax-level3.owl#>
             select ?s where {{ ?s rdf:type biopax:MolecularInteraction . }}
             """)
        return [str(row[0]) for row in qres]



if __name__ == "__main__":
    res = PSICQUICDataSource("P0AEZ3").protein()
    print(res, len(res))
    res = PSICQUICDataSource("P0AEZ3").small_molecule()
    print(res, len(res))
    res = PSICQUICDataSource("P0AEZ3").interaction()
    print(res, len(res))
