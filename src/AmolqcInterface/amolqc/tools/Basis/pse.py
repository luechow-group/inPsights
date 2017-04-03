#
# pse module
#


pse = ['X','H','He',
       'Li','Be','B','C','N','O','F','Ne',
       'Na','Mg','Al','Si','P','S','Cl','Ar',
       'K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr']

pseDict = {'H':'Hydrogen','He':'Helium',
   'Li':'Lithium','Be':'Beryllium','B':'Boron','C':'Carbon','N':'Nitrogen','O':'Oxygen','F':'Fluorine','Ne':'Neon', 
   'Na':'Sodium','Mg':'Magnesium','Al':'Aluminum','Si':'Silicon','P':'Phosphorous','S':'Sulfur','Cl':'Chlorine','Ar':'Argon',
   'K':'Kalium','Ca':'Calcium','Sc':'Scandium','V':'Vanadium','Cr':'Chromium','Mn':'Manganese','Fe':'Iron','Co':'Cobalt','Ni':'Nickel',
   'Cu':'Copper','Zn':'Zinc','Ga':'Gallium','Ge':'Germanium','As':'Arsenic','Se':'Selenium','Br':'Bromine','Kr':'Krypton'}


def elem(idx):
    assert idx>0 and idx<len(pse),"pse.getElement: internal pse too small"
    return pse[idx]

def Z(elem):
    assert elem in pse,"pse.getZ: element is not in internal pse"
    return pse.index(elem)

def elementName(elem):
    assert elem in pseDict.keys(),"pse.getElementName: element not found"
    return pseDict[elem]

if __name__ == "__main__":
    print "testing the pse module"
    assert pse[8] == 'O', "Test 1 failed"
    print Z('Cl')
    print elementName('C')

    
