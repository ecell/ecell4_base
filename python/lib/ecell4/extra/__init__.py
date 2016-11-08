from ecell4.core import Species
try:
    from urllib2 import Request, urlopen
except ImportError:
    from urllib.request import Request, urlopen

class DataSource:

    def __init__(self):
        pass

    @classmethod
    def description(cls, uid):
        return None

class UniProtSource(DataSource):

    def __init__(self):
        pass

    @classmethod
    def getid(cls, obj):
        if isinstance(obj, Species) and obj.has_attribute("uniprot.id"):
            return obj.get_attribute("uniprot.id")
        elif isinstance(obj, str):
            return obj
        else:
            return None

    @classmethod
    def description(cls, obj):
        uid = cls.getid(obj)
        if uid is None:
            return None
        url = 'http://www.uniprot.org/uniprot/{}.txt'.format(uid)
        req = Request(url)
        response = urlopen(req)
        data = response.read()
        return data.decode('utf-8')

def description(obj, database="uniprot"):
    if database == "uniprot":
        return UniProtSource.description(obj)
    return None


if __name__ == "__main__":
    sp = Species("MinD")
    sp.set_attribute("uniprot.id", "P0AEZ3")
    print(description(sp, database="uniprot"))