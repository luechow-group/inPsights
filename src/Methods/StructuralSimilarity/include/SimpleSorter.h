//
// Created by Michael Heuer on 19.06.18.
//

#ifndef AMOLQCPP_SIMPLESORTER_H
#define AMOLQCPP_SIMPLESORTER_H

#include <vector>

#include "StructuralSimilarity.h"

class MolecularSpectrum;

class SimpleSorter{
public:
    std::vector<std::vector<unsigned >> sort(std::vector<MolecularSpectrum> spectra, double threshold = 0.98);
};



#endif //AMOLQCPP_SIMPLESORTER_H
