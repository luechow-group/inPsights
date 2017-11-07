//
// Created by Michael Heuer on 07.11.17.
//

#include <iostream>

#include "RefFileImporter.h"
#include "OptimizationPathFileImporter.h"
#include "ElementInfo.h"


int main(int argc, char *argv[]) {
    std::cout << "Hello" << std::endl;


    RefFileImporter importer("Ethane-max.ref");

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

    //OptimizationPathFileImporter("")

};