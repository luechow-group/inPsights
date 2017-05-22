#include <iostream>
#include "pse.h"
#include "FileXyzInput.h"
#include "Molecule.h"
#include "TestElectronAssigner.h"
#include "HungarianElectronAssigner.h"


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
    input.readElectronCoreAssignations(newMolecule.getCores(),hea);
    input.printAssignations();
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
    input.readElectronCoreAssignations(newMolecule.getCores(),tea);
    input.printAssignations();
    return 0;
}


int testTestElectronAssigner() {
    Molecule newMolecule;
    FileXyzInput input("../testinput/EPA.ref","../testinput/EPA.xyz");
    input.readMoleculeCores(newMolecule);
    TestElectronAssigner tea;
    input.readElectronCoreAssignations(newMolecule.getCores(),tea);
    input.printAssignations();
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
            "LocalSpinMultiplicity reffile xyzfile numOfalphaElectrons bmax/bstart atomnumber StatOnly\n"
            "bmax for assignment based on maximum positions, bstart for assignment based on starting positions\n"
            "atomnumber is the number of the atom for Local Spin Calculation" << std::endl;
}

int main(int argc, char **argv) {
    if(argc!=7){
        show_usage();
        return 1;
    }
    bool basedOnMax=false;
    if(strcmp(argv[4],"bmax")==0)basedOnMax=true;
    bool onlyStat=false;
    if(strcmp(argv[6],"StatOnly")==0)onlyStat=true;
    Molecule newMolecule;
    FileXyzInput input(argv[1],argv[2]);
    input.readMoleculeCores(newMolecule);
    HungarianElectronAssigner hea;
    if(basedOnMax) {
        input.readElectronCoreAssignations(newMolecule.getCores(), hea);
    }
    SpinDeterminer sd(atoi(argv[3]));
    int LocalSpinAtomNumber=atoi(argv[5]);
    if(!onlyStat)std::cout << "SpinQZTotalMolecule\tSpinQZonlyAtomNumber"<<LocalSpinAtomNumber << std::endl;
    std::vector<std::pair<int,int> > SpinQuantumNumbers;
    while(!input.readElectronStructure(newMolecule, sd,basedOnMax?0:&hea)) {
        int TotalSpinQuantumNumber=newMolecule.getTotalSpinQuantumNumber();
        int LocalSpinQuantumNumber=newMolecule.getLocalSpinQuantumNumber(LocalSpinAtomNumber);
        if(!onlyStat) {
            std::cout << TotalSpinQuantumNumber << '\t'
                      << LocalSpinQuantumNumber << std::endl;
        }
        std::vector<std::pair<int,int> >::iterator i;
        for(i=SpinQuantumNumbers.begin();i!=SpinQuantumNumbers.end();i++) {
            if (LocalSpinQuantumNumber == i->first){
                i->second++;
                break;
            }
        }
        if(i==SpinQuantumNumbers.end()) {
            SpinQuantumNumbers.emplace_back(LocalSpinQuantumNumber,1);
        }
    }
    std::cout << "Statistics\n"
            "SpinQuantumNumber Count" << std::endl;
    for(std::vector<std::pair<int,int> >::iterator i=SpinQuantumNumbers.begin();i!=SpinQuantumNumbers.end();i++){
        std::cout << i->first << '\t' << i->second << std::endl;
    }
    std::cout << "Like a Charm!" << std::endl;
    return 0;
}