from bs4 import BeautifulSoup
import requests

def biomodels(biomodels_id):
    xml_soup = BeautifulSoup(requests.get("https://www.ebi.ac.uk/biomodels-main/download?mid=" + biomodels_id).content, 'xml')
    reactions = xml_soup.listOfReactions.find_all('reaction')
    rpks = []
    for r in reactions:
        rpk = []
        reactant = r.listOfReactants.find('speciesReference')['species']
        product = r.listOfProducts.find('speciesReference')['species']
        kinetics = {}
        for k in r.listOfParameters.find_all('parameter'):
            kinetics[k['id']] = k['value']
        rpks.append([reactant, product, kinetics])
    return rpks
