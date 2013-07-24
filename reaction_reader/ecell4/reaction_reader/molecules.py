
import copy
import types
import numbers
import warnings
import functools

#import parseobj

import ecell4.core
#============================================================
# Species Meta Infomation CLASS XXX
class Subunit:
    def __init__(self, name = None):
        self.name = name
        self.modification_list = []
        self.species = None
    def add_modification(self, new_modification):
        self.modification_list.append(new_modification)
    def get_modification(self):
        return self.modification_list

    # XXX Following 2 functions( (set|get)_species ) may not be necessary.
    #   Because it makes cross refenrence between Subunit and Meta_Species, 
    #       the structure of class dependency will be complicated.
    def set_species(self, sp):
        self.species = sp
    def get_species(self):
        return self.species
    def enum_binding_subunit(self):
        ret = []
        for mod in self.modification_list:
            for binding in mod.get_binding():
                ret.append(binding.subunit)
        return ret
    def name(self):
        pass

class Modification:
    def __init__(self, subunit, name, attribute = None):
        self.subunit = subunit  # reference 
        self.name = name
        self.binding = []   # references to binding partner's Modification object.
        self.attribute = attribute 

    def set_binding(self, substrates):
        self.binding = [s for s in substrates if s is not self]
    def get_binding(self):
        return self.binding
    def get_subunit(self):
        return self.subunit
    def name(self):
        return self.name

class Meta_Species(ecell4.core.Species):
    def __init__(self, name):
        self.subunit_list = []
        ecell4.core.Species.__init__(self, name)
    def add_subunit(self, sub):
        self.subunit_list.append(sub)
        sub.set_species(self)
    def get_subunit(self):
        return self.subunit_list

    def dot_output(self):
        acc = []
        print "digraph %s {" % self.name()
        print "\tgraph [label = \"%s\", labelloc = t];" % self.name()
        for sub in self.subunit_list:
            for s in sub.enum_binding_subunit():
                print "\t\"%s\" -> \"%s\";" % (sub.name, s.name)
        print "}"


def enum_edges(meta_species):

    def make_edge(obj1, obj2):
        edge_tuple = ()
        if isinstance(obj1, Subunit):
            if isinstance(obj2, Subunit):
                if obj1.name > obj2.name:
                    edge_tuple = (obj1, obj2)
                else:
                    edge_tuple = (obj2, obj1)
            else:
                edge_tuple = (obj1, obj2)
        else:
            if isinstance(obj2, Subunit):
                edge_tuple = (obj2, obj1)
            else:
                if obj1.name > obj2.name:
                    edge_tuple = (obj1, obj2)
                else:
                    edge_tuple = (obj2, obj1)
        return edge_tuple

    edges = []
    for subunit in meta_species.get_subunit():
        # Subunit - Mofification
        for mod in subunit.modification_list:
            edges.append( make_edge(subunit, mod) )
            # Modification - Modification (binding)
            for bind in mod.get_binding():
                edges.append( make_edge(mod, bind) )
    return edges

def desc(obj):
    s = ""
    if isinstance(obj, Subunit):
        s = "Subunit: " + obj.name
    else:
        s = "Modofication: " + obj.name
    return s

def test_main():
    # mapk(phos=pYT^1).kk(bs^1)
    sp = Meta_Species("mapk")
    subunit_mapk = Subunit("mapk")
    phos_modification = Modification(sp, "YT", "phos")
    subunit_mapk.add_modification(phos_modification)
    sp.add_subunit(subunit_mapk)

    subunit_kk = Subunit("kk")
    bs_modification = Modification(sp, "bs")
    subunit_kk.add_modification(bs_modification)
    sp.add_subunit(subunit_kk)

    bind = [phos_modification, bs_modification]

    for m in bind:
        m.set_binding(bind)
    print "ok"
    #import ipdb; ipdb.set_trace()
    import pdb; pdb.set_trace()
    edges = enum_edges(sp)
    for i in edges:
        print "(", desc(i[0]), " -- ", desc(i[1]), ")"

test_main()
