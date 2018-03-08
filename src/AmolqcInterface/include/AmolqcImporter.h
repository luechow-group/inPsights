//
// Created by Michael Heuer on 06.11.17.
//

#ifndef AMOLQCGUI_AMOLQCIMPORTER_H
#define AMOLQCGUI_AMOLQCIMPORTER_H

#include "Importer.h"

class AmolqcImporter : public Importer{
public:
    explicit AmolqcImporter(const std::string& filename);


    PositionCollection importPositionCollectionBlock(unsigned long startLineIdx,
                                                     unsigned long startLineElement,
                                                     unsigned long numberOfPositions) const;
    
    std::vector<SubstructureDataEntry> countSubstructures(unsigned long startLineIdx,
                                                          unsigned long blockLength) const;

protected:
    SpinTypeCollection getSpinTypeCollection(unsigned long numberOfAlphaElectrons,
                                             unsigned long numberOfBetaElectrons) const;
    
};

#endif //AMOLQCGUI_AMOLQCIMPORTER_H
