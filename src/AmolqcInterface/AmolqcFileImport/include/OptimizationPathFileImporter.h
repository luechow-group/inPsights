// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_OPTIMIZATIONPATHFILEIMPORTER_H
#define INPSIGHTS_OPTIMIZATIONPATHFILEIMPORTER_H

#include "AmolqcImporter.h"
#include <ParticlesVectorCollection.h>

class PathElementDataEntry{
    PathElementDataEntry(unsigned long startLineIdx, unsigned long numberOfPathElements)
            : startLineIdx_(startLineIdx),
              numberOfPathElements_(numberOfPathElements){};
public:
    unsigned long startLineIdx_, numberOfPathElements_;
};

class OptimizationPathFileImporter : public AmolqcImporter{
public:
    explicit OptimizationPathFileImporter(const std::string& filename);

    ElectronsVectorCollection getPath(unsigned long k) const;

    AtomsVector getAtomsVector() const;

    PositionsVector getPositionsVector(unsigned long k, unsigned long m) const;

    unsigned long getNumberOfPaths() const;

private:
    unsigned long calculateLine(unsigned long k, unsigned long m) const;

    AtomsVector atoms_;
    unsigned long numberOfNuclei_,beginOfElectronPositionBlocks_, numberOfElectrons_,
            numberOfAlphaElectrons_,
            numberOfBetaElectrons_,
            numberOfPaths_;
    std::vector<SubstructureDataEntry> substructuresData_;
};

#endif //INPSIGHTS_OPTIMIZATIONPATHFILEIMPORTER_H
