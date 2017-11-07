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

#include "RefFileImporter.h"
#include "OptimizationPathFileImporter.h"

int main(int argc, char const *argv[]) {

    //std::string filename = "Ethane-em-5.wf";
    //WaveFunctionParser waveFunctionParser(filename);
    //waveFunctionParser.readNuclei();
    //auto ac = waveFunctionParser.getAtomCollection();

    std::string filename = "Ethane-max.ref";
    RefFileImporter importer(filename);

    auto ac2 = importer.getAtomCollection();
    for (int i = 0; i < ac2.numberOfParticles(); ++i) {
        std::cout << Elements::ElementInfo::symbol(ac2.elementType(i)) << ac2[i].position().transpose() << std::endl;
    }

    //std::vector<int> a = {};//{6,8,10};
    //a.insert(a.begin()+1,2);
    //std::cout << a[0] << " " << a[1]<< " " << a[2]<< " " << a[3] << std::endl;
    auto ecs = importer.getAllSubstructures(1);
    for (int i = 0; i < ecs.getNumberOfParticleCollections(); ++i) {
        std::cout << ecs.getSpinTypeCollection().spinTypesAsEigenVector().transpose() << std::endl;
        std::cout << ecs[i].positionsAsEigenVector().transpose() << std::endl;
    }

    std::string filename2 = "Diborane-Paths.300";
    OptimizationPathFileImporter importer2(filename2,1);

    auto ecs2 = importer2.getPath(1);
    for (int i = 0; i < ecs2.getNumberOfParticleCollections(); ++i) {
        std::cout << ecs2.getSpinTypeCollection().spinTypesAsEigenVector().transpose() << std::endl;
        std::cout << ecs2[i].positionsAsEigenVector().transpose() << std::endl;
    }

}
