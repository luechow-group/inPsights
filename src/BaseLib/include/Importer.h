//
// Created by Michael Heuer on 02.11.17.
//

#ifndef AMOLQCGUI_IMPORTER_H
#define AMOLQCGUI_IMPORTER_H

#include <fstream>
#include <tuple>
#include "AtomCollection.h"
#include "ElectronCollections.h"


class Importer{

public:
    explicit Importer(const std::string& filename);
    ~Importer();

    std::string getLine(unsigned long idx) const;

    template<typename Out>
    void split(const std::string &s, char delimiter, Out result) const;


    std::vector<std::string> split(const std::string &s) const;
    std::vector<std::string> split(const std::string &s, char delim) const;



private:
    template<class OutIt>
    void read_lines(std::istream& is, OutIt dest);

    std::string filename_;
    std::ifstream file_;
    std::vector<std::string> lines_;
};



class Format{

public:

};

class RefFileImporter : public Importer{
public:
    RefFileImporter(const std::string& filename);

    AtomCollection getAtomCollection();

    ElectronCollection getElectronCollection(unsigned long k, unsigned long m);
    unsigned long getNumberOfMaxima(unsigned long k, unsigned long m);


private:

    void countSubstructures();
    unsigned long calculateLine(unsigned long k, unsigned long m);

    unsigned long numberOfNuclei_,numberOfElectrons_, numberOfAlphaElectrons_,
            numberOfSuperstructures_, totalNumberOfMaxima_, maximalNumberOfSubstructures;

    // line idx, numerOfSubstructures, totalNumberOfMaxima
    std::vector<std::tuple<unsigned long,unsigned long,unsigned long>> substructuresData_;
};

class WfFileImporter : public Importer{
public:
    WfFileImporter(const std::string& filename);


private:

};

#endif //AMOLQCGUI_IMPORTER_H
