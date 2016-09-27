import re

from rdflib import Namespace
from rdflib.namespace import RDF, RDFS, SKOS
from rdflib.term import URIRef

try:
    from . import rdf
    from . import pdb
except SystemError:
    import rdf
    import pdb

def description(entity):
    desc = []

    entity_id = UniProtDataSource.parse_entity(entity)
    if entity_id is not None:
        src = UniProtDataSource(entity_id)
        if src.obsolete():
            desc.append("UniProtKB - {} (This entry is obsolete)".format(entity_id))
            desc.append("URL: {}".format(UniProtDataSource.link(entity)))
            return desc
        desc.append("UniProtKB - {} ({})".format(entity_id, ', '.join(src.mnemonic())))
        desc.append("Protein: {}".format(', '.join(src.structured_name())))
        desc.append("Gene: {}".format(', '.join(src.gene())))
        desc.append("Organism: {}".format(', '.join(src.organism())))
        desc.append("Function: {}".format(', '.join(src.function_annotation())))
        desc.append("URL: {}".format(UniProtDataSource.link(entity)))
        return desc

    entity_id = UniProtLocationDataSource.parse_entity(entity)
    if entity_id is not None:
        src = UniProtLocationDataSource(entity_id)
        desc.append("UniProtKB - {} ({})".format(', '.join(src.pref_label()), src.get_type()))
        desc.append("Definition: {}".format(', '.join(src.comment())))
        desc.append("Synonyms: {}".format(', '.join(src.alt_label())))
        desc.append("PartOf: {}".format(', '.join(src.part_of())))
        desc.append("GO: {}".format(', '.join(src.go())))
        desc.append("URL: {}".format(UniProtLocationDataSource.link(entity)))
        return desc

    return desc

def whereis(entity):
    desc = []
    entity_id = UniProtDataSource.parse_entity(entity)
    if entity_id is not None:
        src = UniProtDataSource(entity_id)
        for component, topology in src.subcellular_location():
            tmp = description(component)
            if len(tmp) > 0 and len(desc) > 0:
                desc.append('')
            desc.extend(tmp)
            if topology is None:
                continue
            tmp = description(topology)
            if len(tmp) > 0 and len(desc) > 0:
                desc.append('')
            desc.extend(tmp)
    return desc

class UniProtDataSourceBase(rdf.RDFDataSourceBase):

    GRAPH = {}
    UNIPROT = Namespace("http://purl.uniprot.org/core/")
    UPDB = Namespace("http://purl.uniprot.org/database/")

    def __init__(self, url=None, cache=True):
        rdf.RDFDataSourceBase.__init__(self, url, cache)

class UniProtTaxonDataSource(UniProtDataSourceBase):

    URL = "http://www.uniprot.org/taxonomy/{entry_id}.rdf"

    def __init__(self, entry_id=None, url=None, cache=True):
        if url is not None:
            UniProtDataSourceBase.__init__(
                self, url, cache)
        elif entry_id is not None:
            UniProtDataSourceBase.__init__(
                self, self.URL.format(entry_id=entry_id), cache)
        else:
            UniProtDataSourceBase.__init__(self, None, cache)

    def scientific_name(self):
        return [str(obj) for obj in self.graph.objects(predicate=self.UNIPROT.scientificName)]

class UniProtLocationDataSource(UniProtDataSourceBase):

    URL = "http://www.uniprot.org/locations/{entity_id}.rdf"

    def __init__(self, entity=None, url=None, cache=True):
        if url is not None:
            UniProtDataSourceBase.__init__(
                self, url, cache)
        elif entity is not None:
            entity_id = self.parse_entity(entity)
            if entity_id is not None:
                UniProtDataSourceBase.__init__(
                    self, self.URL.format(entity_id=entity_id), cache)
            else:
                UniProtDataSourceBase.__init__(self, None, cache)
        else:
            UniProtDataSourceBase.__init__(self, None, cache)

    @classmethod
    def parse_entity(cls, entity):
        idpttrn= r'(?P<prefix>SL-)?[0-9]{1,4}'
        uri1 = r'http://purl.uniprot\.org/locations/(?P<id>{})(\.rdf)?'.format(idpttrn)
        uri2 = r'http://www.uniprot\.org/locations/(?P<id>{})(\.rdf)?'.format(idpttrn)

        if isinstance(entity, str):
            mobj = re.match(r'^{}$'.format(idpttrn), entity)
            if mobj is not None:
                if mobj.group('prefix') is None:
                    return 'SL-{:04d}'.format(int(entity))
                else:
                    return entity
            mobj = re.match(uri1, entity)
            if mobj is not None:
                if mobj.group('prefix') is None:
                    return 'SL-{:04d}'.format(int(mobj.group('id')))
                else:
                    return mobj.group('id')
            mobj = re.match(uri2, entity)
            if mobj is not None:
                if mobj.group('prefix') is None:
                    return 'SL-{:04d}'.format(int(mobj.group('id')))
                else:
                    return mobj.group('id')
        # else:
        #     import ecell4
        #     if isinstance(entity, ecell4.Species) and entity.has_attribute(collection):
        #         return cls.parse_entity(entity.get_attribute(collection))
        return None  #XXX: Error

    @classmethod
    def link(cls, entity):
        entity_id = cls.parse_entity(entity)
        assert entity_id is not None
        return 'http://www.uniprot.org/locations/{}'.format(entity_id)

    def pref_label(self):
        return [str(obj) for obj in self.graph.objects(predicate=SKOS.prefLabel)]

    def alt_label(self):
        return [str(obj) for obj in self.graph.objects(predicate=SKOS.altLabel)]

    def alias(self):
        return [str(obj) for obj in self.graph.objects(predicate=self.UNIPROT.alias)]

    def comment(self):
        return [str(obj) for obj in self.graph.objects(predicate=RDFS.comment)]

    def go(self):
        return [str(sub) for sub in self.graph.subjects(predicate=self.UNIPROT.database, object=self.UPDB.go)]

    def part_of(self):
        return [str(obj) for obj in self.graph.objects(predicate=self.UNIPROT.partOf)]

    def get_type(self):
        qres = self.query(
            """prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            prefix skos: <http://www.w3.org/2004/02/skos/core#>
            select ?type where
            {{
            ?s
            rdf:type ?type ;
            skos:prefLabel ?label .
            }}
            """)
        uri = [str(row[0]) for row in qres]
        assert len(uri) == 1
        label = [str(obj) for obj in rdf.RDFDataSourceBase(uri[0]).graph.objects(subject=URIRef(uri[0]), predicate=RDFS.label)]
        assert len(label) == 1
        return label[0]

class UniProtDataSource(UniProtDataSourceBase):

    URL = "http://www.uniprot.org/uniprot/{entity_id}.rdf"

    def __init__(self, entity=None, url=None, cache=True):
        if url is not None:
            UniProtDataSourceBase.__init__(
                self, url, cache)
        elif entity is not None:
            entity_id = self.parse_entity(entity)
            if entity_id is not None:
                UniProtDataSourceBase.__init__(
                    self, self.URL.format(entity_id=entity_id), cache)
            else:
                UniProtDataSourceBase.__init__(self, None, cache)
        else:
            UniProtDataSourceBase.__init__(self, None, cache)

    @classmethod
    def parse_entity(cls, entity):
        # http://www.ebi.ac.uk/miriam/main/collections/MIR:00000005
        collection = 'uniprot'
        idpttrn = r'([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9])(\.\d+)?'
        uri1 = r'http://identifiers\.org/uniprot/(?P<id>{})'.format(idpttrn)
        uri2 = r'http://www.uniprot\.org/uniprot/(?P<id>{}).rdf'.format(idpttrn)
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
        return 'http://www.uniprot.org/uniprot/{}'.format(entity_id)

    def mnemonic(self):
        return [str(obj) for obj in self.graph.objects(predicate=self.UNIPROT.mnemonic)]  #XXX: Specify its subject

    def obsolete(self):
        return any(self.graph.objects(predicate=self.UNIPROT.obsolete))  #XXX: Specify its subject

    def gene(self):
        return [str(obj) for obj in self.objects(self.UNIPROT.Gene, SKOS.prefLabel)]

    def locus_name(self):
        return [str(obj) for obj in self.objects(self.UNIPROT.Gene, self.UNIPROT.locusName)]

    def function_annotation(self):
        return [str(obj) for obj in self.objects(self.UNIPROT.Function_Annotation, RDFS.comment)]

    def simple_sequence(self):
        return [str(obj) for obj in self.objects(self.UNIPROT.Simple_Sequence, RDF.value)]

    def sequence_annotation(self, uri=UniProtDataSourceBase.UNIPROT.Sequence_Annotation):
        # http://www.uniprot.org/core/
        # http://www.uniprot.org/help/sequence_annotation
        qres = rdf.RDFDataSourceBase(str(self.UNIPROT)).query(
            """prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>
            prefix uniprot: <http://purl.uniprot.org/core/>
            select ?s ?label ?comment ?see_also where
            {{
            {{
            {{
            ?s
            rdfs:subClassOf+ <{uri}> ;
            rdfs:label ?label .
            }}
            union
            {{
            ?s rdfs:label ?label .
            filter( ?s = <{uri}> ).
            }}
            optional {{ ?s rdfs:comment ?comment }}
            optional {{ ?s rdfs:seeAlso ?see_also . }}
            }}
            }}
            """.format(uri=str(uri)))
        names = {}
        for row in qres:
            name, label = str(row[0]), str(row[1])
            value = dict(name=name, label=label)
            if row[2] is not None:
                value['comment'] = str(row[2])
            if row[3] is not None:
                value['see_also'] = str(row[3])
            names[name] = value

        qres = self.query("""prefix rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>
            prefix faldo: <http://biohackathon.org/resource/faldo#>
            prefix uniprot: <http://purl.uniprot.org/core/>
            select ?type ?begin ?end ?comment where
            {{
            ?s
            rdf:type ?type ;
            uniprot:range ?range .
            optional {{ ?s rdfs:comment ?comment }} .
            filter( ?type in ({}) ) .
            ?range faldo:begin ?begin_ ; faldo:end ?end_ .
            ?begin_ faldo:position ?begin .
            ?end_ faldo:position ?end .
            }}
            """.format(', '.join('<{}>'.format(name) for name in names.keys())))
        retval = []
        for row in qres:
            name, begin, end = str(row[0]), int(row[1]), int(row[2])
            value = dict(begin=begin, end=end, type=names[name])
            if row[3] is not None:
                value['comment'] = str(row[3])
            retval.append(value)
        return retval

    def molecule_processing(self):
        return self.sequence_annotation(uri=self.UNIPROT.Molecule_Processing_Annotation)

    def region(self):
        return self.sequence_annotation(uri=self.UNIPROT.Region_Annotation)

    def site(self):
        return self.sequence_annotation(uri=self.UNIPROT.Site_Annotation)

    def modification(self):
        return self.sequence_annotation(uri=self.UNIPROT.Modification_Annotation)

    def natural_variation(self):
        return self.sequence_annotation(uri=self.UNIPROT.Natural_Variation_Annotation)

    def experimental_information(self):
        return self.sequence_annotation(uri=self.UNIPROT.Experimental_Information_Annotation)

    def secondary_structure(self):
        return self.sequence_annotation(uri=self.UNIPROT.Secondary_Structure_Annotation)

    def subcellular_location(self):
        qres = self.query("""
            prefix uniprot: <http://purl.uniprot.org/core/>
            select ?cellular_component ?topology where
            {{
            ?s
            rdf:type uniprot:Subcellular_Location_Annotation ;
            uniprot:locatedIn ?w .
            ?w
            uniprot:cellularComponent ?cellular_component .
            optional {{ ?w uniprot:topology ?topology }} .
            }}
            """)
        retval = []
        for row in qres:
            if row[1] is None:
                retval.append((str(row[0]), None))
            else:
                retval.append((str(row[0]), str(row[1])))
        return retval

    def organism(self):
        retval = []
        for obj1 in self.graph.objects(predicate=self.UNIPROT.organism):
            mobj = re.match("http:\/\/purl\.uniprot\.org\/taxonomy\/([0-9]+)", str(obj1))
            if mobj is None:
                continue
            # taxon_id = mobj.group(1)
            # retval.extend(UniProtTaxonDataSource(taxon_id).scientific_name())
            retval.extend(UniProtTaxonDataSource(url=str(obj1)).scientific_name())
        return retval

    def structured_name(self):
        return [str(obj) for obj in self.objects(self.UNIPROT.Structured_Name, self.UNIPROT.fullName)]

    def structure_resource(self):
        return [str(sub) for sub in self.graph.subjects(predicate=RDF.type, object=self.UNIPROT.Structure_Resource)]

    def pdb(self):
        retval = []
        for sub in self.graph.subjects(predicate=RDF.type, object=self.UNIPROT.Structure_Resource):
            if URIRef("http://purl.uniprot.org/database/PDB") not in self.graph.objects(subject=sub, predicate=self.UNIPROT.database):
                continue
            # mobj = re.match("http:\/\/rdf\.wwpdb\.org\/pdb\/([0-9A-Za-z]+)", str(sub))
            # if mobj is None:
            #     continue
            # pdb_id = mobj.group(1).upper()
            retval.extend(pdb.PDBDataSource(url=str(sub)).identifier())
        return retval

    def database(self, dbname):
        return [str(sub) for sub in self.graph.subjects(predicate=self.UNIPROT.database, object=self.UPDB[dbname])]

    def biogrid(self):
        return self.database("BioGrid")


if __name__ == "__main__":
    # print(description("P0AEZ3"))
    print(description("http://identifiers.org/uniprot/P0AEZ3"))
    # print(description("http://www.uniprot.org/uniprot/P0AEZ3.rdf"))

    print(UniProtDataSource("P0AEZ3").gene()[0])
    print(UniProtDataSource("P0AEZ3").locus_name())
    print(UniProtDataSource("P0AEZ3").function_annotation()[0])
    print(UniProtDataSource("P0AEZ3").organism()[0])
    print(UniProtDataSource("P0AEZ3").structure_resource())
    print(UniProtDataSource("P0AEZ3").pdb())
    print(UniProtDataSource("P0AEZ3").biogrid())
    print(UniProtDataSource("P0AEZ3").database("IntAct"))

    print(UniProtDataSource("P28482").molecule_processing())
    print(UniProtDataSource("P28482").region())
    print(UniProtDataSource("P28482").site())
    print(UniProtDataSource("P28482").modification())
    print(UniProtDataSource("P28482").natural_variation())
    print(UniProtDataSource("P28482").experimental_information())
    print(UniProtDataSource("P28482").secondary_structure())

    print(UniProtDataSource("P28482").subcellular_location())
    print(whereis("P28482"))
