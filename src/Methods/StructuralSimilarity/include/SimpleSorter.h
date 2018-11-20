//
// Created by Michael Heuer on 19.06.18.
//

#ifndef INPSIGHTS_SIMPLESORTER_H
#define INPSIGHTS_SIMPLESORTER_H

#include <vector>

#include "StructuralSimilarity.h"

class MolecularSpectrum;

class SimpleSorter{
public:
    std::vector<std::vector<unsigned >> sort(std::vector<MolecularSpectrum> spectra, double threshold = 0.98);
};



#endif //INPSIGHTS_SIMPLESORTER_H
