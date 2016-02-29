from libsbml import *

reader = SBMLReader()

class BioModels(object):
    """Class to represent a BioModels model."""
    def __init__(self, modelid):
        self.modelid = modelid
        document = reader.readSBML(modelid)
        self.model = document.getModel()
        self.reactions = [self.model.getReaction(x) for x in range(0, self.model.getNumReactions())]
