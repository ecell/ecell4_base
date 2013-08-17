import species
import string   # for convert into .bngl

def convert2bng_species(self):
    return ".".join([subunit.convert2bng() for subunit in self.subunits])

def convert2bng_subunit(self):
    mods1, mods2 = [], []
    for mod, (state, binding) in self.modifications.items():
        if state == '':
            if binding == '_':
                mods1.append("%s!+" % (mod))
            elif binding != '':
                mods1.append("%s!%s" % (mod, binding))
            else:
                mods1.append(mod)
        else:
            if binding == '_':
                mods2.append("%s~%s!+" % (mod, binding))
            elif binding == "":
                mods2.append("%s~%s" % (mod, state))
            else:
                mods2.append("%s~%s!%s" % (mod, state, binding))
    mods1.sort()
    mods2.sort()
    mods1.extend(mods2)
    #return str(self).translate(string.maketrans('=^', '~!'))
    return "%s(%s)" % (self.name, ",".join(mods1))

def convert2bng_reactionrule(self):
    return "%s -> %s" % (
        "+".join([sp.convert2bng() for sp in self.reactants()]),
        "+".join([sp.convert2bng() for sp in self.products()]))

def convert2bng_include_reactants(self):
    return "inlude_reactants(%d,%s)" % (self.__idx, self.__pttrn)
def convert2bng_exclude_reactants(self):
    return "exlude_reactants(%d,%s)" % (self.__idx, self.__pttrn)
def convert2bng_include_products(self):
    return "inlude_products(%d,%s)" % (self.__idx, self.__pttrn)

def convert2bng_exclude_products(self):
    return "exlude_products(%d,%s)" % (self.__idx, self.__pttrn)


class Convert2BNGManager(object):
    def __init__(self, species, rules):
        self.__species = species
        self.__rules = rules
        self.__modification_collection_dict = {}
        if 0 < len(species) and 0 < len(rules):
            self.build_modification_collection_dict()

        self.initialize_methods()

    def initialize_methods(self):
        species.Species.convert2bng = convert2bng_species
        species.Subunit.convert2bng = convert2bng_subunit
        species.ReactionRule.convert2bng = convert2bng_reactionrule
        species.IncludeReactants.convert2bng = convert2bng_include_reactants
        species.ExcludeReactants.convert2bng = convert2bng_exclude_reactants
        species.IncludeProducts.convert2bng = convert2bng_include_products
        species.ExcludeProducts.convert2bng = convert2bng_exclude_products 


    def write_section_seed_species(self, fd):
        fd.write("begin seed species\n")
        for i, (sp, attr) in enumerate( self.__species ):
            fd.write("\t%s\t%f\n" % (sp.convert2bng(), attr))
        fd.write("end seed species\n")

    def write_section_molecule_types(self, fd):
        def build_molecules_type_query_list(current_dict):
            retval = []
            for su_name in current_dict:
                mod_list = []
                for m, state_list in current_dict[su_name].items():
                    mod = "%s" % m
                    for state in list(set(state_list)):
                        if state != '' and state != "_":
                            mod = "%s~%s" % (mod, state)
                    mod_list.append(mod)
                retval.append("%s(%s)" % (su_name, ','.join(mod_list) ))
            return retval

        # write
        fd.write("begin molecule types\n")
        for s in build_molecules_type_query_list(self.__modification_collection_dict):
            fd.write("\t%s\n" % s)
        fd.write("end molecule types\n")

    def write_section_reaction_rules(self, fd):
        fd.write("begin reaction rules\n")
        for i, rr in enumerate(self.__rules):
            s = "\t%s\t%f" % (rr.convert2bng(),rr.options()[0])
            #fd.write("\t%s\t%f" % (rr.convert2bng(),rr.options()[0]))
            for cond in rr.options():
                if isinstance(cond, species.Option):
                    s = "%s %s" % (s, cond.convert2bng() )
            s += "\n"
            fd.write(s)
        fd.write("end reaction rules\n")

    def build_modification_collection_dict(self):
        def add_modification_collection_dict_subunit(current_dict, subunit):
            su_name = subunit.get_name()
            if not current_dict.has_key( su_name ):
                current_dict[subunit.get_name()] = {}
            for mod, (state, binding) in subunit.get_modifications_list().items():
                if mod in current_dict[su_name]:
                    current_dict[su_name][mod].append(state)
                else:
                    current_dict[su_name][mod] = [state]
            return current_dict
        temp_dict = {}
        reactants = []
        products = []
        for rr in self.__rules:
            reacntants = rr.reactants()
            products = rr.products()
            for r in reactants:
                for su in r.get_subunit_list():
                    temp_dict = add_modification_collection_dict_subunit(temp_dict, su)
            for p in products:
                for su in p.get_subunit_list():
                    temp_dict = add_modification_collection_dict_subunit(temp_dict, su)
        self.__modification_collection_dict = temp_dict
    
    def get_modification_collection_dict(self):
        return self.__modification_collection_dict
