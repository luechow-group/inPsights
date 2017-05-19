//
// Created by Moria on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_ELECTRONPERMUTATION_H
#define LOCALSPINMULTIPLICITY_ELECTRONPERMUTATION_H
#include "Electron.h"
#include <vector>


class ElectronPermutation {
public:
    ElectronPermutation(std::vector<int> permutation);
    virtual ~ElectronPermutation();
    std::vector<Electron>& permutate(std::vector<Electron>& electrons);
private:
    std::vector<int> permutation;
};


#endif //LOCALSPINMULTIPLICITY_ELECTRONPERMUTATION_H
