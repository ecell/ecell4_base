import re

try:
    from urllib.request import Request, urlopen, HTTPError  # Python3
except ImportError:
    from urllib2 import Request, urlopen, HTTPError  # Python2

from xml.dom import minidom


def read_url(url):
    f = urlopen(url)
    content = f.read().decode('utf-8')
    f.close()
    try:
        f = urlopen(url)
        content = f.read().decode('utf-8')
        f.close()
    except IOError:
        #XXX: say something here
        content = None
    return content

def description(entity):
    entity_id = PubMedDataSource.parse_entity(entity)
    if entity_id is not None:
        entry = []
        src = PubMedDataSource(entity)
        entry.append(('PubMed', '{}'.format(entity_id), ' - '))
        entry.append(('Title', src.data['Title']))
        entry.append(('Author(s)', ', '.join(src.data['AuthorList'])))
        entry.append(('Source', src.data['Source']))
        entry.append(('SO', src.data['SO']))
        entry.append(('URL', src.link(entity_id)))
        return [entry]

    return []

class PubMedDataSource(object):

    def __init__(self, entity=None):
        self.entity = entity

        if entity is not None:
            entity_id = self.parse_entity(entity)
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={}".format(entity_id)
            data = self.parse_esummary(read_url(url))
            assert len(data) == 1
            self.data = data[0]
        else:
            self.data = None

    @classmethod
    def parse_entity(cls, entity):
        # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000015
        collection = 'pubmed'
        idpttrn = r'\d+'
        uri1 = r'https://www.ncbi.nlm.nih.gov/pubmed/(?P<id>{})'.format(idpttrn)
        uri2 = r'http://identifiers.org/pubmed/(?P<id>{})'.format(idpttrn)
        if isinstance(entity, str):
            if re.match(r'^{}$'.format(idpttrn), entity) is not None:
                return entity
            mobj = re.match(uri1, entity)
            if mobj is not None:
                return mobj.group('id')
            mobj = re.match(uri2, entity)
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
        return 'https://www.ncbi.nlm.nih.gov/pubmed/{}'.format(entity_id)

    @classmethod
    def parse_esummary(cls, esummary):
        retval = []
        doc = minidom.parseString(esummary)
        for entry_node in doc.getElementsByTagName('DocSum'):
            entry = {}
            entry['ID'] = entry_node.getElementsByTagName('Id')[0].firstChild.data
            for item in entry_node.getElementsByTagName('Item'):
                name = item.getAttribute('Name')
                if name in ('Title', 'Volume', 'Issue', 'Pages', 'Source', 'PubDate', 'SO'):
                    entry[name] = item.firstChild.data
                elif name == 'AuthorList':
                    entry['AuthorList'] = [author.firstChild.data for author in item.getElementsByTagName('Item') if author.getAttribute('Name') == 'Author']
            retval.append(entry)
        return retval


if __name__ == "__main__":
    print(description("8752322"))
