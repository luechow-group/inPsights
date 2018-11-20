//
// Created by Michael Heuer on 06.11.17.
//

#ifndef INPSIGHTS_AMOLQCIMPORTER_H
#define INPSIGHTS_AMOLQCIMPORTER_H

#include "Importer.h"
#include <PositionsVector.h>
#include <TypesVector.h>

class AmolqcImporter : public Importer{
public:
    explicit AmolqcImporter(const std::string& filename);


    PositionsVector importPositionsVectorBlock(unsigned long startLineIdx,
                                               unsigned long startLineElement,
                                               unsigned long numberOfPositions) const;
    
    std::vector<SubstructureDataEntry> countSubstructures(unsigned long startLineIdx,
                                                          unsigned long blockLength) const;

protected:
    SpinTypesVector getSpinTypesVector(unsigned long numberOfAlphaElectrons,
                                       unsigned long numberOfBetaElectrons) const;
    
};

#endif //INPSIGHTS_AMOLQCIMPORTER_H
