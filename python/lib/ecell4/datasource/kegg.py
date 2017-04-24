import requests
from Bio.KEGG import Compound
import io

def description(entity):
    if entity is not None:
        entry = []
        kegg_get = requests.get('http://rest.kegg.jp/get/' + entity.get_attribute('kegg'))
        records = Compound.parse(io.StringIO(kegg_get.content.decode("utf-8")))
        #print(len(list(records)))
        record = list(records)[0]
        #print(record.name)
        entry.append(('KEGG_NAMES', ', '.join(record.name)))
        return [entry]

