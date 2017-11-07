//
// Created by Michael Heuer on 06.11.17.
//

#ifndef AMOLQCGUI_AMOLQCIMPORTER_H
#define AMOLQCGUI_AMOLQCIMPORTER_H

#include "Importer.h"



class AmolqcImporter : public Importer{
public:
    AmolqcImporter(const std::string& filename);


    ParticleCollection importParticleCollectionBlock(unsigned long startLineIdx,
                                                     unsigned long startLineElement,
                                                     unsigned long numberOfParticles) const;
    SpinTypeCollection getSpinTypeCollection(unsigned long numberOfAlphaElectrons,
                                             unsigned long numberOfBetaElectrons) const;

    std::vector<SubstructureDataEntry> countSubstructures(unsigned long startLineIdx,
                                                          unsigned long blockLength) const;



};

#endif //AMOLQCGUI_AMOLQCIMPORTER_H
