//
// Created by Michael Heuer on 26.01.17.
//

#include <iostream>
#include "ElectronicWaveFunction.h"
#include "ElementInfo.h"

int main(int argc, char const *argv[]) {

    auto &wf = ElectronicWaveFunction::getInstance("Ethane-em-5.wf");

    auto atomCollection = wf.getAtomCollection();

    std::cout << atomCollection << std::endl;

    auto atom1 = atomCollection.atom(1);

    std::cout << atom1 << std::endl
              << "symbol = "<< atom1.elementType() << std::endl
              << "mass   = "<< Elements::ElementInfo::mass(atom1.elementType()) << std::endl
              << std::endl;

    Elements::ElementType fluorineElementType = Elements::ElementType::F;
    std::cout << "symbol = "<< fluorineElementType << std::endl
              << "p-electrons = "<< Elements::ElementInfo::pElectrons(fluorineElementType)
              << std::endl;
}
