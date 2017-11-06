//
// Created by Michael Heuer on 09.01.17.
//

#ifndef MAXREFPARSER_H
#define MAXREFPARSER_H

#include <fstream>

#include "AtomCollection.h"

class WaveFunctionParser {
public:
    WaveFunctionParser(const std::string &filename);
    void readNuclei();

    ~WaveFunctionParser();

    AtomCollection getAtomCollection() { return atomCollection_; };

private:
    std::string filename_;
    std::ifstream file_;

    AtomCollection atomCollection_;
};

#endif //MAXREFPARSER_H
