// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_REFFILEIMPORTER_H
#define INPSIGHTS_REFFILEIMPORTER_H

#include "AmolqcImporter.h"
#include <ParticlesVectorCollection.h>

class RefFileImporter : public AmolqcImporter{
public:
    RefFileImporter(const std::string& filename);

    AtomsVector getAtomsVector();

    SpinTypesVector getSpinTypesVector() const;
    PositionsVector getPositionsVector(unsigned long k, unsigned long m) const;
    ElectronsVector getMaximaStructure(unsigned long k, unsigned long m) const;
    ElectronsVectorCollection getAllSubstructures(unsigned long k) const;

    unsigned long getNumberOfMaxima(unsigned long k, unsigned long m) const;
    double getNegativeLogarithmizedProbabilityDensity(unsigned long k, unsigned long m) const;

    unsigned long numberOfSuperstructures();

private:
    unsigned long calculateLine(unsigned long k, unsigned long m) const;

    unsigned long numberOfNuclei_,
            numberOfElectrons_,
            numberOfAlphaElectrons_,
            numberOfBetaElectrons_,
            numberOfSuperstructures_,
            totalNumberOfMaxima_,
            maximalNumberOfSubstructures; ;
    const unsigned numberOfLinesAboveCoordinatesBlock = 2;

    std::vector<SubstructureDataEntry> substructuresData_;
};

#endif //INPSIGHTS_REFFILEIMPORTER_H
