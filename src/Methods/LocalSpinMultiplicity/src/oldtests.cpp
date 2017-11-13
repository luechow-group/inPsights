//
// Created by Moria on 20.07.2017.
//

/*

#include <iostream>
#include "pse.h"
#include "FileXyzInput.h"
#include "Molecule.h"
#include "TestElectronAssigner.h"
#include "HungarianElectronAssigner.h"
#include "SpinQuantumNumberCounter.h"

int testReadElectronStructureOwnAssignment() {
    const std::vector<int> CH1={0,1};
    const std::vector<int> CH2={0,2};
    const std::vector<int> CH3={0,3};
    const std::vector<int> CH4={0,4};
    Molecule newMolecule;
    FileXyzInput input("../resources/EPA.ref","../resources/EPA.xyz");
    input.readMoleculeCores(newMolecule);
    HungarianElectronAssigner hea;
    SpinDeterminer sd(5);
    std::cerr << "SpinQZTotalMolecule\tSpinQZonlyC\tSpinQZfragmentCandH1\tSpinQZfragmentCandH2\tSpinQZfragmentCandH3\tSpinQZfragmentCandH4\t" << std::endl;
    while(!input.readElectronStructure(newMolecule, sd, &hea)){
        std::cerr << newMolecule.getTotalSpinQuantumNumber() << '\t'
                  << newMolecule.getLocalSpinQuantumNumber(0) << '\t'
                  << newMolecule.getLocalSpinQuantumNumber(CH1) << '\t'
                  << newMolecule.getLocalSpinQuantumNumber(CH2) << '\t'
                  << newMolecule.getLocalSpinQuantumNumber(CH3) << '\t'
                  << newMolecule.getLocalSpinQuantumNumber(CH4) << std::endl;
    }
    return 0;
}

int testReadElectronStructure() {
    const std::vector<int> CH1={0,1};
    const std::vector<int> CH2={0,2};
    const std::vector<int> CH3={0,3};
    const std::vector<int> CH4={0,4};
    Molecule newMolecule;
    FileXyzInput input("../resources/EPA.ref","../resources/EPA.xyz");
    input.readMoleculeCores(newMolecule);
    HungarianElectronAssigner hea;
    input.readElectronCoreAssignments(newMolecule.getCores(),hea);
    input.printAssignments();
    SpinDeterminer sd(5);
    std::cerr << "SpinQZTotalMolecule\tSpinQZonlyC\tSpinQZfragmentCandH1\tSpinQZfragmentCandH2\tSpinQZfragmentCandH3\tSpinQZfragmentCandH4\t" << std::endl;
    while(!input.readElectronStructure(newMolecule, sd)) {
        std::cerr << newMolecule.getTotalSpinQuantumNumber() << '\t'
                  << newMolecule.getLocalSpinQuantumNumber(0) << '\t'
                  << newMolecule.getLocalSpinQuantumNumber(CH1) << '\t'
                  << newMolecule.getLocalSpinQuantumNumber(CH2) << '\t'
                  << newMolecule.getLocalSpinQuantumNumber(CH3) << '\t'
                  << newMolecule.getLocalSpinQuantumNumber(CH4) << std::endl;
    }
    return 0;
}

int testHungarianElectronAssigner() {
    Molecule newMolecule;
    FileXyzInput input("../resources/EPA.ref","../resources/EPA.xyz");
    input.readMoleculeCores(newMolecule);
    HungarianElectronAssigner tea;
    input.readElectronCoreAssignments(newMolecule.getCores(),tea);
    input.printAssignments();
    return 0;
}


int testTestElectronAssigner() {
    Molecule newMolecule;
    FileXyzInput input("../resources/EPA.ref","../resources/EPA.xyz");
    input.readMoleculeCores(newMolecule);
    TestElectronAssigner tea;
    input.readElectronCoreAssignments(newMolecule.getCores(),tea);
    input.printAssignments();
    return 0;
}

int testFileXyzInput(){
    Molecule newMolecule;
    FileXyzInput input("../resources/EPA.ref","../resources/EPA.xyz");
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

 */
