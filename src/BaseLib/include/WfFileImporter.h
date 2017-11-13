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
    WfFileImporter(const std::string &filename);
    void readNuclei();


    std::pair<bool, unsigned long> findTag(const std::string &tag,
                                           unsigned long startLine = 0);




        AtomCollection getAtomCollection(); //{ return atomCollection_; };

private:
    void readGeneralBlock();

    bool angstromQ_;
    int charge_;
    unsigned long numberOfNuclei_,multiplicity_;
    std::string basis_, title_, jastrow_;
};

#endif //WFFILEIMPORTER_H
