import species
import string   # for convert into .bngl
import copy

# Label Related
def is_label(s):
    return 1 < len(s) and s[0] == '_'

def check_label(self):
    for su in self.subunits:
        for mod, (state, binding) in su.get_modifications_list().items():
            pass
        pass

def check_labeled_subunit(su):
    retval = []
    for mod, (state, binding) in su.get_modifications_list().items():
        if 1 < len(state) and state[0] == '_':
            retval.append( (state,mod) )    # state should be label.
    return retval

def check_label_containing_reaction(rr):
    def add_label_dict(current_dict, label, subunit, mod):
        if label not in current_dict:
            current_dict[label] = [(subunit.get_name(), mod)]
        else:
            current_dict[label].append( (subunit.get_name(), mod) )
        return current_dict

    reactant_labels = {}
    product_labels = {}
    for r in rr.reactants():
        for su in r.get_subunit_list():
            d = check_labeled_subunit(su)
            for (label, mod) in d:
                reactant_labels = add_label_dict(reactant_labels, label, su, mod)
    for p in rr.products():
        for su in p.get_subunit_list():
            d = check_labeled_subunit(su)
            for (lebel, mod) in d:
                product_labels = add_label_dict(product_labels, label, su, mod)
    return (reactant_labels, product_labels)

def build_label_expanded_reactionrule(rr):
    (reactant_labels, product_labels) = check_label_containing_reaction(rr)
    print "reactant labels"
    print reactant_labels
    print "product_labels"
    print product_labels
        
# Formatters for each class
# class Species
def convert2bng_species(self, labels = None):
    return ".".join([subunit.convert2bng(labels) for subunit in self.subunits])

# class Subunit
def convert2bng_subunit(self, labels = None):
    mods1, mods2 = [], []
    for mod, (state, binding) in self.modifications.items():
        if state == '':
            if binding == '_':
                mods1.append("%s!+" % (mod))
            elif binding != '':
                mods1.append("%s!%s" % (mod, binding))
            else:
                mods1.append(mod)
        elif is_label(state):
            if labels != None and (state in labels):
                if binding == '_':
                    mods2.append("%s~%s!+" % (mod, binding))
                elif binding == "":
                    mods2.append("%s~%s" % (mod, labels[state]))
                else:
                    mods2.append("%s~%s!%s" % (mod, labels[state], binding))
            else:
                print ("Warning: candidates for label %s was not found" % state)
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

# class ReactionRule
def convert2bng_reactionrule(self, labels = None):
    return "%s -> %s" % (
        "+".join([sp.convert2bng(labels) for sp in self.reactants()]),
        "+".join([sp.convert2bng(labels) for sp in self.products()]))

# classe Options
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
        for rr in self.__rules:
            rr_expanded = self.build_label_expanded_reactionrule(rr)
            if isinstance(rr_expanded, list):
                for rr_query in rr_expanded:
                    s = "\t%s\t%f" % (rr_query, rr.options()[0])
                    for cond in rr.options():
                        if isinstance(cond, species.Option):
                            s = "%s %s" % (s, cond.convert2bng() )
                    s += "\n"
                    fd.write(s)
            else:
                s = "\t%s\t%f" % (rr_expanded, rr.options()[0])
                for cond in rr.options():
                    s = "%s %s" % (s, cond.convert2bng() )
                s += "\n"
                fd.write(s)

        '''
        for i, rr in enumerate(self.__rules):
            s = "\t%s\t%f" % (rr.convert2bng(),rr.options()[0])
            #fd.write("\t%s\t%f" % (rr.convert2bng(),rr.options()[0]))
            for cond in rr.options():
                if isinstance(cond, species.Option):
                    s = "%s %s" % (s, cond.convert2bng() )
            s += "\n"
        '''
        fd.write("end reaction rules\n")

    def build_modification_collection_dict(self):
        def add_modification_collection_dict_subunit(current_dict, subunit):
            su_name = subunit.get_name()
            if not current_dict.has_key( su_name ):
                current_dict[subunit.get_name()] = {}
            for mod, (state, binding) in subunit.get_modifications_list().items():
                if not is_label(state):
                    if mod in current_dict[su_name]:
                        current_dict[su_name][mod].append(state)
                    else:
                        current_dict[su_name][mod] = [state]
            return current_dict

        temp_dict = {}
        # Build modification dictionary by species
        for (sp, attr) in self.__species:
            for su in sp.get_subunit_list():
                temp_dict = add_modification_collection_dict_subunit(temp_dict, su)
        # Build modification dictionary by ReactionRules
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

    def build_label_expanded_reactionrule(self, rr):
        def convert2bngl_label_expanded_reactionrule_recursive(
                rr, label_list, candidates, acc, bnglstr_acc):
            acc_conbination = copy.deepcopy(acc)
            for mod in candidates[ label_list[0] ]:
                acc_conbination[label_list[0]] = mod
                if len(label_list) == 1:
                    bnglstr_acc.append( rr.convert2bng(acc_conbination) )
                else:
                    convert2bngl_label_expanded_reactionrule_recursive(
                            rr, label_list[1:], candidates, acc_conbination, bnglstr_acc)

        (reactant_labels, product_labels) = check_label_containing_reaction(rr)
        if reactant_labels or product_labels:
            bngl_strs = []
            modification_candidates = {}
            for label, pos in reactant_labels.items():
                modification_candidates[label] = []
                for (su, mod) in pos:
                    modification_candidates[label].extend(
                            self.__modification_collection_dict[su][mod])
                if 0 < len(modification_candidates):
                    print ("%s candidates" % (label))
                    print modification_candidates[label]
        
            convert2bngl_label_expanded_reactionrule_recursive(
                    rr, modification_candidates.keys(), modification_candidates, {}, bngl_strs)
            print "lebel expanded_rules"
            for s in bngl_strs:
                print s
            return bngl_strs
        else: #containing no labels.
            return rr.convert2bng()

