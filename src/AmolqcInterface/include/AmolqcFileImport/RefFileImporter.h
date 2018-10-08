//
// Created by Michael Heuer on 06.11.17.
//

#ifndef AMOLQCPP_REFFILEIMPORTER_H
#define AMOLQCPP_REFFILEIMPORTER_H

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

#endif //AMOLQCPP_REFFILEIMPORTER_H
