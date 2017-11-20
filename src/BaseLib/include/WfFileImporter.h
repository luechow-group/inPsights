//
// Created by Michael Heuer on 09.01.17.
//

#ifndef WFFILEIMPORTER_H
#define WFFILEIMPORTER_H

#include <fstream>
#include "AmolqcImporter.h"
#include "AtomCollection.h"

enum JastrowTypes{
    none,
    IC,
    SM
};

class WfFileImporter : public AmolqcImporter {
public:
    explicit WfFileImporter(const std::string &filename);

    std::pair<bool, unsigned long> findTag(const std::string &tag,
                                           unsigned long startLine = 0);

    AtomCollection getAtomCollection();
    SpinTypeCollection getSpinTypeCollection();
    unsigned long getNumberOfElectrons();
    unsigned long getNumberOfAlphaElectrons();
    unsigned long getNumberOfBetaElectrons();
    long getCharge(){ return charge_; };
    unsigned long getMultiplicity(){ return multiplicity_; };
    std::string getBasisSet(){ return basis_; };

private:
    void readGeneralBlock();

    bool bohrQ_;
    unsigned long charge_, numberOfNuclei_,multiplicity_;
    std::string basis_, title_, jastrow_;
};

#endif //WFFILEIMPORTER_H
