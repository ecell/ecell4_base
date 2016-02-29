from libsbml import *

reader = SBMLReader()

class BioModels(object):
    """Class to represent a BioModels model."""
    def __init__(self, modelid):
        self.modelid = modelid
        document = reader.readSBML(modelid)
        self.model = document.getModel()
