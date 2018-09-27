//
// Created by Michael Heuer on 07.11.17.
//

#ifndef AMOLQCPP_OPTIMIZATIONPATHFILEIMPORTER_H
#define AMOLQCPP_OPTIMIZATIONPATHFILEIMPORTER_H

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
    OptimizationPathFileImporter(const std::string& filename, unsigned long  multiplicity);

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

#endif //AMOLQCPP_OPTIMIZATIONPATHFILEIMPORTER_H
