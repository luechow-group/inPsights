//
// Created by Michael Heuer on 26.01.17.
//

#include <iostream>
#include <iomanip>

#include "ChemicalSystem.h"
#include "WfFileImporter.h"
#include <ElementInfo.h>

#include "RefFileImporter.h"
#include "OptimizationPathFileImporter.h"

#include "CollectionParser.h"

int main(int argc, char const *argv[]) {

    std::string filename = "Ethane-em-5.wf";
    WfFileImporter waveFunctionParser(filename);
    auto ac = waveFunctionParser.getAtomCollection();

    filename = "Ethane-max.ref";
    RefFileImporter importer(filename);

    auto ac2 = importer.getAtomCollection();
    for (int i = 0; i < ac2.numberOfParticles(); ++i) {
        std::cout << Elements::ElementInfo::symbol(ac2.elementType(i)) << ac2[i].position().transpose() << std::endl;
    }

    auto ecs = importer.getAllSubstructures(1);
    for (int i = 0; i < ecs.length(); ++i) {
        std::cout << ecs.getSpinTypeCollection().spinTypesAsEigenVector().transpose() << std::endl;
        std::cout << ecs[i].positionsAsEigenVector().transpose() << std::endl;
    }

    std::string filename2 = "Diborane-Paths.300";
    OptimizationPathFileImporter importer2(filename2,1);

    auto ecs2 = importer2.getPath(1);
    for (int i = 0; i < ecs2.length(); ++i) {
        std::cout << ecs2.getSpinTypeCollection().spinTypesAsEigenVector().transpose() << std::endl;
        std::cout << ecs2[i].positionsAsEigenVector().transpose() << std::endl;
    }


    RefFileImporter refFileImporter("Ethane-max.ref");
    auto ec = refFileImporter.getMaximaStructure(1,1);
    CollectionParser collectionParser;

    auto acJson = collectionParser.atomCollectionToJson(ac);
    collectionParser.writeJSON(acJson,"ac.json");

    auto ecJson = collectionParser.electronCollectionToJson(ec);
    collectionParser.writeJSON(ecJson,"ec.json");

    auto acp = collectionParser.atomCollectionFromJson("ac.json");
    auto ecp = collectionParser.electronCollectionFromJson("ec.json");

    std::cout << acp.positionsAsEigenVector().transpose() << std::endl;
    std::cout << static_cast<ElementTypeCollection>(acp).elementTypesAsEigenVector().transpose() << std::endl;

    std::cout << ecp.positionsAsEigenVector().transpose() << std::endl;
    std::cout << static_cast<SpinTypeCollection>(ecp).spinTypesAsEigenVector().transpose() << std::endl;

}
