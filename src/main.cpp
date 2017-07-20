#include <iostream>
#include "pse.h"
#include "FileXyzInput.h"
#include "Molecule.h"
#include "TestElectronAssigner.h"
#include "HungarianElectronAssigner.h"
#include "SpinQuantumNumberCounter.h"


/*int testReadElectronStructureOwnAssignment() {
    const std::vector<int> CH1={0,1};
    const std::vector<int> CH2={0,2};
    const std::vector<int> CH3={0,3};
    const std::vector<int> CH4={0,4};
    Molecule newMolecule;
    FileXyzInput input("../testinput/EPA.ref","../testinput/EPA.xyz");
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
    FileXyzInput input("../testinput/EPA.ref","../testinput/EPA.xyz");
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
    FileXyzInput input("../testinput/EPA.ref","../testinput/EPA.xyz");
    input.readMoleculeCores(newMolecule);
    HungarianElectronAssigner tea;
    input.readElectronCoreAssignments(newMolecule.getCores(),tea);
    input.printAssignments();
    return 0;
}


int testTestElectronAssigner() {
    Molecule newMolecule;
    FileXyzInput input("../testinput/EPA.ref","../testinput/EPA.xyz");
    input.readMoleculeCores(newMolecule);
    TestElectronAssigner tea;
    input.readElectronCoreAssignments(newMolecule.getCores(),tea);
    input.printAssignments();
    return 0;
}

int testFileXyzInput(){
    Molecule newMolecule;
    FileXyzInput input("../testinput/EPA.ref","../testinput/EPA.xyz");
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
}*/

void show_usage()
{
    std::cout << "Usage: \n"
            "LocalSpinMultiplicity reffile xyzfile numOfalphaElectrons bmax/bstart StatOnly atomnumbers \n"
            "bmax for assignment based on maximum positions, bstart for assignment based on starting positions\n"
            "atomnumbers are the numbers of the atoms for Local Spin Calculation\n Given none, whole molecule will be looked at." << std::endl;
}

int main(int argc, char **argv) {
    if(argc<6){
        show_usage();
        return 1;
    }

    bool basedOnMax=false;
    if(strcmp(argv[4],"bmax")==0)basedOnMax=true;

    bool onlyStat=false;
    if(strcmp(argv[5],"StatOnly")==0)onlyStat=true;

    std::vector<int> fragmentAtomNumbers;
    for(int i=6;i<argc;i++){
        fragmentAtomNumbers.emplace_back(atoi(argv[i]));
    }

    Molecule newMolecule;
    FileXyzInput input(argv[1],argv[2]);
    input.readMoleculeCores(newMolecule);
    HungarianElectronAssigner hea;
    if(basedOnMax) {
        input.readElectronCoreAssignments(newMolecule.getCores(), hea);
    }
    SpinDeterminer sd(atoi(argv[3]));
    std::cout << (basedOnMax?"bmax ":"bstart ");
    switch(fragmentAtomNumbers.size()){
        case 0:
            std::cout << "SpinQZwholeMolecule" << std::endl;
            break;
        case 1:
            std::cout << "SpinQZonlyAtomNumber "<< fragmentAtomNumbers[0] << std::endl;
            break;
        default:
            std::cout << "SpinQZonlyAtomNumbers ";
            for(std::vector<int>::iterator i=fragmentAtomNumbers.begin();i!=fragmentAtomNumbers.end();i++){
                std::cout << *i << ' ';
            }
            std::cout << std::endl;
            break;
    }
    SpinQuantumNumberCounter SQNCounter;
    while(!input.readElectronStructure(newMolecule, sd,basedOnMax?0:&hea)) {
        if(!onlyStat) {
            std::cout << static_cast<double>(SQNCounter.addNumber(newMolecule.getLocalSpinQuantumNumber(fragmentAtomNumbers)))/2 << std::endl;
        } else {
            SQNCounter.addNumber(newMolecule.getLocalSpinQuantumNumber(fragmentAtomNumbers));
        }
    }
    SQNCounter.printStatsSpinQuantumNumber();
    SQNCounter.printStatsMultiplicities();
    std::cout << "Like a Charm!" << std::endl;
    return 0;
}