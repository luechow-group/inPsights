//
// Created by Michael Heuer on 06.11.17.
//

#ifndef AMOLQCGUI_REFFILEIMPORTER_H
#define AMOLQCGUI_REFFILEIMPORTER_H

#include "AmolqcImporter.h"

class RefFileImporter : public AmolqcImporter{
public:
    RefFileImporter(const std::string& filename);

    AtomCollection getAtomCollection();

    SpinTypesVector getSpinTypesVector() const;
    PositionsVector getPositionsVector(unsigned long k, unsigned long m) const;
    ElectronCollection getMaximaStructure(unsigned long k, unsigned long m) const;
    ElectronsVectorCollection getAllSubstructures(unsigned long k) const;

    unsigned long getNumberOfMaxima(unsigned long k, unsigned long m) const;
    double getNegativeLogarithmizedProbabilityDensity(unsigned long k, unsigned long m) const;

private:
    unsigned long calculateLine(unsigned long k, unsigned long m) const;

    unsigned long numberOfNuclei_,
            numberOfElectrons_,
            numberOfAlphaElectrons_,
            numberOfBetaElectrons_,
            numberOfSuperstructures_, totalNumberOfMaxima_, maximalNumberOfSubstructures;
    // line idx, numerOfSubstructures, totalNumberOfMaxima

    //std::vector<std::tuple<unsigned long,unsigned long,unsigned long>> substructuresData_;
    std::vector<SubstructureDataEntry> substructuresData_;
};

#endif //AMOLQCGUI_REFFILEIMPORTER_H
