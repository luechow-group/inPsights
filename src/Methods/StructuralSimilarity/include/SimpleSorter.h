// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_SIMPLESORTER_H
#define INPSIGHTS_SIMPLESORTER_H

#include <vector>
#include "StructuralSimilarity.h"
#include "MolecularSpectrum.h"

class SimpleSorter{
public:
    std::vector<std::vector<unsigned >> sort(std::vector<SOAP::MolecularSpectrum> spectra, double threshold = 0.98);
};



#endif //INPSIGHTS_SIMPLESORTER_H
