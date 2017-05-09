# coding: utf-8

__author__ = "Keita Sasaki"

import json

try:
    from urllib2 import Request, urlopen
except ImportError:
    from urllib.request import Request, urlopen


class biogridDataSource(object):

    histry = {}

    def __init__(self, ACCESS_KEY):
        self.ak = ACCESS_KEY

    def query(self, ACCESS_POINT, geneList=None, organisms=None, cache=True):
        geneList = geneList or []
        organisms = organisms or []

        options = ""
        if geneList:
            options += "&geneList={}".format("|".join(geneList))
        if organisms:
            options += "&taxid={}".format("|".join(organisms))

        url = "http://webservice.thebiogrid.org/{AP}/?format=json&searchNames=true{OPTIONS}&accesskey={AK}".format(AP=ACCESS_POINT, OPTIONS=options, AK=self.ak)

        if not cache or url not in self.histry.keys():
            req = Request(url)
            response = urlopen(req)
            data = response.read().decode("utf-8")

            if cache:
                self.histry[url] = data

        else:
            data = self.histry[url]

        return json.loads(data)

    def organisms(self):
        return self.query("organisms").values()

    def orgmaker(self, org=None):
        org = org or []

        organisms_reverse = dict(
            [(value, key) for key, value in self.query("organisms").items()])
        if org:
            for i in range(len(org)):
                if isinstance(org[i], int):
                    org[i] = str(org[i])
                elif isinstance(org[i], str):
                    try:
                        org[i] = organisms_reverse[org[i]]
                    except KeyError:
                        pass
        return org

    def interactions(self, geneList=None, org=None):
        geneList = geneList or []
        org = org or []

        organisms = self.query("organisms")
        org = self.orgmaker(org)

        querydata = self.query("interactions", geneList, org).values()
        returnData = []
        for i in querydata:
            if i["OFFICIAL_SYMBOL_A"] in geneList:
                dataDict = {"symA": {"name": i["OFFICIAL_SYMBOL_A"],
                                     "biogridID": i["BIOGRID_ID_A"],
                                     "organism": {"name": organisms[str(i["ORGANISM_A"])], "id": i["ORGANISM_A"]}},
                            "symB": {"name": i["OFFICIAL_SYMBOL_B"],
                                     "biogridID": i["BIOGRID_ID_B"],
                                     "organism": {"name": organisms[str(i["ORGANISM_B"])], "id": i["ORGANISM_B"]}},
                            "pubmedID": i["PUBMED_ID"]
                            }
            else:
                dataDict = {"symA": {"name": i["OFFICIAL_SYMBOL_B"],
                                     "biogridID": i["BIOGRID_ID_B"],
                                     "organism": {"name": organisms[str(i["ORGANISM_B"])], "id": i["ORGANISM_B"]}},
                            "symB": {"name": i["OFFICIAL_SYMBOL_A"],
                                     "biogridID": i["BIOGRID_ID_A"],
                                     "organism": {"name": organisms[str(i["ORGANISM_A"])], "id": i["ORGANISM_A"]}},
                            "pubmedID": i["PUBMED_ID"]
                            }
            returnData.append(dataDict)

        return returnData

    def interactor(self, geneList=None, org=None):
        """
        Supposing geneList returns an unique item.
        """
        geneList = geneList or []
        organisms = organisms or []

        querydata = self.interactions(geneList, org)
        returnData = {}
        for i in querydata:
            if not returnData.get(i["symB"]["name"]):
                returnData[i["symB"]["name"]] = {"interactions": []}
            returnData[i["symB"]["name"]]["interactions"].append(i)
        return returnData


if __name__ == "__main__":
    ACCESSKEY = "YOUR_ACCESSKEY"
    interactions = biogridDataSource(ACCESSKEY).interactions(["MDM2"], [9606])
    interactor = biogridDataSource(ACCESSKEY).interactor(["MDM2"])
    print(interactor.keys())
