#include <iostream>
#include "pse.h"
#include "FileXyzInput.h"
#include "Molecule.h"

int testFileXyzInput(){
    Molecule newMolecule;
    FileXyzInput input("../input/EPA.ref","../input/EPA.xyz");
    input.readMoleculeCores(newMolecule);
    return 0;
}

int testFindElement() {
    std::string testElement[11] = {"C", "H", "Xe", "Pb", "Kr", "Ni", "Cl", "At", "Ra", "Lu", "Hg"};
    int testElementOT[11] = {6,1,54,82,36,28,17,85,88,71,80};
    for (int i = 0; i < 11; i++) {
        int OZ=Pse::findElement(testElement[i]);
        std::cout << testElement[i] << ": " << OZ << " should be: " << testElementOT[i] << std::endl;
        if(OZ!=testElementOT[i]){
            std::cout << "Test failed" << std::endl;
            return 1;
        }

    }
    return 0;
}

int main() {
    if(testFindElement())return 1;
    if(testFileXyzInput())return 1;
    return 0;
}