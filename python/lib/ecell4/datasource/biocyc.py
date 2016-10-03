import re

try:
    from urllib.request import Request, urlopen, HTTPError  # Python3
except ImportError:
    from urllib2 import Request, urlopen, HTTPError  # Python2

from xml.dom import minidom
import xml.dom


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

class BioCycDataSource(object):

    def __init__(self, entity=None):
        if entity is not None:
            entity_id = self.parse_entity(entity)
            orgid, frameid = entity_id.split(':')
            url = "http://websvc.biocyc.org/getxml?id={0}:{1}&detail=low".format(orgid, frameid)
            data = self.parse_ptools_xml(read_url(url))
            assert len(data) == 1
            self.data = data
        else:
            self.data = None

    @classmethod
    def parse_entity(cls, entity):
        collection = 'biocyc'
        idpttrn = r'[A-Z-0-9]+(?<!CHEBI)(\:)?[A-Za-z0-9+_.%-]+'
        uri1 = r'http://identifiers.org/biocyc/(?P<id>{})'.format(idpttrn)
        if isinstance(entity, str):
            if re.match(r'^{}$'.format(idpttrn), entity) is not None:
                return entity
            mobj = re.match(uri1, entity)
            if mobj is not None:
                return mobj.group('id')
            # mobj = re.match(uri2, entity)
            # if mobj is not None:
            #     return mobj.group('id')
        else:
            import ecell4
            if isinstance(entity, ecell4.Species) and entity.has_attribute(collection):
                return cls.parse_entity(entity.get_attribute(collection))
        return None  #XXX: Error

    @classmethod
    def link(cls, entity):
        entity_id = cls.parse_entity(entity)
        assert entity_id is not None
        return ""

    @classmethod
    def __parse_ptools_xml(cls, node, tags=None, ignores=None, unique=False):
        entries = []
        for entry_node in node.childNodes:
            entry = {}
            if entry_node.nodeType != xml.dom.Node.ELEMENT_NODE:
                continue
            elif entry_node.tagName == 'coefficient':
                entries[-1]['coefficient'] = entry_node.firstChild.data
                continue
            elif tags is None or entry_node.tagName in tags:
                entry['type'] = entry_node.tagName
                entry['frameid'] = entry_node.getAttribute('frameid')
                entry['orgid'] = entry_node.getAttribute('orgid')

                if entry_node.hasAttribute("class"):
                    value = entry_node.getAttribute("class")
                    if value == 'true':
                        entry["class"] = True
                    elif value == 'false':
                        entry["class"] = False
                    else:
                        raise ValueError("Unknown value for 'class' attribute was given [{}]".format(value))

                for item in entry_node.childNodes:
                    if item.nodeType != xml.dom.Node.ELEMENT_NODE:
                        continue
                    elif item.tagName in ('common-name', ):
                        # item.getAttribute('datatype') == 'string':
                        assert item.tagName not in entry.keys()
                        value = item.firstChild.data.strip()
                        entry[item.tagName] = value
                    elif item.tagName in ('synonym', 'ec-number'):
                        # item.getAttribute('datatype') == 'string':
                        value = item.firstChild.data.strip()
                        if item.tagName not in entry.keys():
                            entry[item.tagName] = [value]
                        else:
                            entry[item.tagName].append(value)
                    elif item.tagName in ('component', 'component-of', 'parent', 'gene', 'left', 'right', 'enzyme', 'reaction', 'instance'):
                        if item.tagName not in entry.keys():
                            entry[item.tagName] = cls.__parse_ptools_xml(item)
                        else:
                            entry[item.tagName].extend(cls.__parse_ptools_xml(item))
                    elif item.tagName in ('catalyzes', 'enzymatic-reaction'):
                        if item.tagName not in entry.keys():
                            entry[item.tagName] = cls.__parse_ptools_xml(item, 'Enzymatic-Reaction')
                        else:
                            entry[item.tagName].extend(cls.__parse_ptools_xml(item, 'Enzymatic-Reaction'))
                    elif item.tagName == 'cml':
                        for subitem in item.childNodes:
                            if subitem.nodeType != xml.dom.Node.ELEMENT_NODE:
                                continue
                            elif subitem.tagName == 'molecule':
                                for subsubitem in subitem.childNodes:
                                    if subsubitem.nodeType != xml.dom.Node.ELEMENT_NODE:
                                        continue
                                    elif subsubitem.tagName == 'formula':
                                        assert 'formula' not in entry.keys()
                                        entry['formula'] = subsubitem.getAttribute('concise')
                                    elif subsubitem.tagName == 'float':
                                        if 'float' not in entry.keys():
                                            entry['float'] = [{'title': subsubitem.getAttribute('title'), 'units': subsubitem.getAttribute('units'), 'value': subsubitem.firstChild.data}]
                                        else:
                                            entry['float'].append({'title': subsubitem.getAttribute('title'), 'units': subsubitem.getAttribute('units'), 'value': subsubitem.firstChild.data})
                                    elif subsubitem.tagName == 'string':
                                        if 'string' not in entry.keys():
                                            entry['string'] = [{'title': subsubitem.getAttribute('title'), 'value': subsubitem.firstChild.data}]
                                        else:
                                            entry['string'].append({'title': subsubitem.getAttribute('title'), 'value': subsubitem.firstChild.data})
            elif ignores is not None and entry_node.tagName in ignores:
                continue
            else:
                raise ValueError('Unknown tag name was given [{}].'.format(entry_node.tagName))
            entries.append(entry)
        if unique:
            assert len(entries) == 1
            return entries[0]
        return entries

    @classmethod
    def parse_ptools_xml(cls, content):
        doc = minidom.parseString(content)
        return cls.__parse_ptools_xml(doc.firstChild, ('Protein', 'Compound', 'Reaction', 'Gene'), ('metadata', ))


if __name__ == "__main__":
    src = BioCycDataSource('ECOLI:GLUCOSE-1-PHOSPHAT-CPLX')
    print(src.data)
    src = BioCycDataSource('ECOLI:EG10597-MONOMER')
    print(src.data)
    src = BioCycDataSource('ECOLI:GLC-1-P')
    print(src.data)
    src = BioCycDataSource('ECOLI:Glucopyranose')
    print(src.data)
    src = BioCycDataSource('ECOLI:GLUCOSE-1-PHOSPHAT-RXN')
    print(src.data)

    src = BioCycDataSource('ECOCYC:EG10597-MONOMER')  #XXX: "EcoCyc:EG10597-MONOMER" doesn't work.
    print(src.data)
    src = BioCycDataSource('ECOL316407:JW1164-MONOMER')
    print(src.data)
