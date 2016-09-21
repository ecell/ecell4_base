try:
    from urllib.request import Request, urlopen, HTTPError  # Python3
except ImportError:
    from urllib2 import Request, urlopen, HTTPError  # Python2

from rdflib.namespace import RDF
from rdflib import Namespace

import xml.dom.minidom

try:
    from . import rdf
except SystemError:
    import rdf



def read_url(url):
    try:
        f = urlopen(url)
        content = f.read()
        f.close()
    except IOError:
        content = None
    return content

def get_active_services():
    url = 'http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=ACTIVE&format=xml'
    content = read_url(url)
    dom = xml.dom.minidom.parseString(content)
    services = []
    for elem in dom.getElementsByTagName('service'):
        name = elem.getElementsByTagName('name')[0].firstChild.data
        restUrl = elem.getElementsByTagName('restUrl')[0].firstChild.data
        services.append((name, restUrl))
    return services

class PSICQUICDataSource(rdf.RDFDataSourceBase):

    BIOPAX = Namespace("http://www.biopax.org/release/biopax-level3.owl#")
    URL = "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/{entry_id}?format=rdf-xml"

    ACTIVE_SERVICES = dict(get_active_services())

    def __init__(self, entry_id=None, cache=True, services=None):
        rdf.RDFDataSourceBase.__init__(self, None, cache)
        self.entry_id = entry_id

        if services is None:
            self.services = tuple(self.ACTIVE_SERVICES.keys())
        elif isinstance(services, str):
            if services in self.ACTIVE_SERVICES.keys():
                self.services = [services]
            else:
                self.services = []
        else:
            self.services = [
                name for name in services if name in self.ACTIVE_SERVICES.keys()]

    def set_graph(self, service_name):
        if self.entry_id is None:
            return

        self.url = "{:s}interactor/{:s}?format=rdf-xml".format(self.ACTIVE_SERVICES[service_name], self.entry_id)

        if self.cache and self.url in self.GRAPH.keys():
            self.graph = self.fetch(self.url, self.cache)
            return

        cnt = int(read_url("{:s}interactor/{:s}?format=count".format(self.ACTIVE_SERVICES[service_name], self.entry_id)))
        # print(service_name, self.url, cnt)

        if cnt == 0:
            if self.cache:
                self.GRAPH[self.url] = None
            self.graph = None
        else:
            try:
                self.graph = self.fetch(self.url, self.cache)
            except HTTPError as e:
                if e.code in (500, 406, 400, 200):
                    if e.code == 500:
                        msg = "HTTP Error {:d}: Internal server error".format(e.code)
                    elif e.code == 406:
                        msg = "HTTP Error {:d}: Format not supported".format(e.code)
                    elif e.code == 400:
                        msg = "HTTP Error {:d}: Too many results for exporting in XML, Biopax and RDF".format(e.code)
                    elif e.code == 200:
                        msg = "HTTP Error {:d}: Not an error. Request is OK".format(e.code)
                    else:
                        msg = e.reason()
                    if self.cache:
                        self.GRAPH[self.url] = None
                    self.graph = None
                    # print(msg)
                else:
                    raise e

    def subjects(self, key):
        retval = []
        for name in self.services:
            self.set_graph(name)
            if self.graph is None:
                continue
            retval.extend(
                [str(sub) for sub in self.graph.subjects(RDF.type, self.BIOPAX[key])])
        return retval

    def protein(self):
        return tuple(set(self.subjects("ProteinReference")))

    def small_molecule(self):
        return tuple(set(self.subjects("SmallMoleculeReference")))

    def interactor(self):
        return tuple(set(self.protein() + self.small_molecule()))

    def interaction(self):
        return tuple(set(self.subjects("MolecularInteraction")))


if __name__ == "__main__":
    print(get_active_services())

    # services = None
    services = "IntAct"
    res = PSICQUICDataSource("P0AEZ3", services=services).protein()
    print(res, len(res))
    res = PSICQUICDataSource("P0AEZ3", services=services).small_molecule()
    print(res, len(res))
    res = PSICQUICDataSource("P0AEZ3", services=services).interaction()
    print(res, len(res))
