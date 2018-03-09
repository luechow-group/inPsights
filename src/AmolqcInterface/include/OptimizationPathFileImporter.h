//
// Created by Michael Heuer on 07.11.17.
//

#ifndef AMOLQCGUI_OPTIMIZATIONPATHFILEIMPORTER_H
#define AMOLQCGUI_OPTIMIZATIONPATHFILEIMPORTER_H

#include "AmolqcImporter.h"


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

    PositionCollection getPositionCollection(unsigned long k, unsigned long m) const;

    unsigned long getNumberOfPaths() const;

private:
    unsigned long calculateLine(unsigned long k, unsigned long m) const;

    unsigned long numberOfElectrons_,
            numberOfAlphaElectrons_,
            numberOfBetaElectrons_,
            numberOfPaths_;
    std::vector<SubstructureDataEntry> substructuresData_;
};

#endif //AMOLQCGUI_OPTIMIZATIONPATHFILEIMPORTER_H
