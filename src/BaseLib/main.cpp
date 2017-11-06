//
// Created by Michael Heuer on 26.01.17.
//

#include <iostream>
#include <iomanip>
#include <string>

#include "ChemicalSystem.h"
#include "WaveFunctionParser.h"
#include "Importer.h"
#include <ElementInfo.h>

int main(int argc, char const *argv[]) {

    std::string filename = "Ethane-em-5.wf";

    WaveFunctionParser waveFunctionParser(filename);
    waveFunctionParser.readNuclei();

    auto ac = waveFunctionParser.getAtomCollection();

    filename = "Ethane-max.ref";

    //Importer importer(filename);
    //std::cout << importer.getLine(10) << std::endl;
    RefFileImporter importer(filename);

    auto ac2 = importer.getAtomCollection();
    for (int i = 0; i < ac2.numberOfParticles(); ++i) {
        std::cout << Elements::ElementInfo::symbol(ac2.elementType(i)) << ac2[i].position().transpose() << std::endl;
    }

}
