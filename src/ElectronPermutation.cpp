//
// Created by Moria on 19.05.2017.
//

#include "ElectronPermutation.h"

ElectronPermutation::ElectronPermutation(std::vector<int> permutation): permutation(permutation){

}

ElectronPermutation::~ElectronPermutation() {

}

std::vector<Electron> &ElectronPermutation::permutate(std::vector<Electron> &electrons) {
    std::vector<Electron> helpElectrons;
    for(int i=0;i<electrons.size();i++){
        helpElectrons.push_back(electrons[permutation[i]]);
    }
    electrons=helpElectrons;
    return electrons;
}
