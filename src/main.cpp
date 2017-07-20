//
// Created by Morian Sonnet on 16.05.2017.
//

#include <iostream>
#include "Molecule.h"
#include "FileXyzInput.h"
#include "HungarianElectronAssigner.h"
#include "SpinProjectionQuantumNumberCounter.h"

void show_usage();

int main(int argc, char **argv) {
    if(argc<6){        //assume user is not aware of usage
        show_usage();
        return 1;       //exit with error code 1
    }

    bool basedOnMax=false; //check for bmode
    if(strcmp(argv[4],"bmax")==0)basedOnMax=true;

    bool onlyStat=false;   //check for StatMode
    if(strcmp(argv[5],"StatOnly")==0)onlyStat=true;

    std::vector<int> fragmentAtomNumbers;   //add atoms taken into account for local spin calculation
    for(int i=6;i<argc;i++){
        fragmentAtomNumbers.emplace_back(atoi(argv[i]));
    }

    Molecule newMolecule;
    FileXyzInput input(argv[1],argv[2]);
    input.readMoleculeCores(newMolecule);
    HungarianElectronAssigner hea;
    if(basedOnMax) {        //in this case the Assignments need to be done only once
        input.readElectronCoreAssignments(newMolecule.getCores(), hea);
    }
    SpinDeterminer sd(atoi(argv[3]));       //pass number of electrons with alpha-spin to SpinDeterminer
    std::cout << (basedOnMax?"bmax ":"bstart ");        //show bmode
    switch(fragmentAtomNumbers.size()){                 //show used Fragment
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
    SpinProjectionQuantumNumberCounter SQNCounter;
    while(!input.readElectronStructure(newMolecule, sd, basedOnMax?nullptr:&hea)) { //for every Electron Arrangement
        if(!onlyStat) {
            std::cout <<
                static_cast<double>(SQNCounter.addNumber(newMolecule.getLocalSpinQuantumNumber(fragmentAtomNumbers)))/2
                      << std::endl;
        } else {
            SQNCounter.addNumber(newMolecule.getLocalSpinQuantumNumber(fragmentAtomNumbers));//add MS value to Statistics
        }
    }
    SQNCounter.printStatsSpinProjectionQuantumNumber();
    SQNCounter.printStatsMultiplicity();
    std::cout << "Like a Charm!" << std::endl; //A wonderful phrase one can grep on
    return 0;
}

void show_usage()
{
    std::cout << "Usage: \n"
        "LocalSpinMultiplicity reffile xyzfile numAlpha bmode StatMode atomnumbers \n"
        "\n"
        "reffile:    Reference File obtained from Electron-Pair-Analysis-Calculation\n"
        "xyzfile:    xyz-File obtained from Electron-Pair-Analysis-Calculation with write_xyz option\n"
        "numAlpha:   Number of Electrons with alpha-Spin\n"
        "bmode:      Electronarrangement the Calculation should be based on\n"
        "               bmax for Maximumarrangements\n"
        "               bstart for Arrangements according to amount square of wavefunction\n"
        "StatMode:   StatOnly for only Statistics of S and Ms. Anything else for MS of every Arrangement\n"
        "atomnumbers:Atoms of Fragment to calculate Local Spin. Given none, whole molecule is taken"
              << std::endl;

}
