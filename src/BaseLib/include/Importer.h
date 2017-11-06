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

class SubstructureDataEntry{
public:
    SubstructureDataEntry(unsigned long startingLine,
                          unsigned long numberOfSubstructures,
                          unsigned long totalNumberOfMaxima)
            : startingLine_(startingLine),
              numberOfSubstructures_(numberOfSubstructures),
              totalNumberOfMaxima_(totalNumberOfMaxima) {};

    unsigned long startingLine_, numberOfSubstructures_, totalNumberOfMaxima_;
};



class AmolqcImporter : public Importer{
public:
    AmolqcImporter(const std::string& filename);


    ParticleCollection importParticleCollectionBlock(unsigned long startLine,
                                                     unsigned long startLineElement,
                                                     unsigned long numberOfParticles) const;
    SpinTypeCollection createSpinTypeCollection(unsigned long numberOfAlphaElectrons,
                                                unsigned long numberOfBetaElectrons) const;
private:
};


class RefFileImporter : public AmolqcImporter{
public:
    RefFileImporter(const std::string& filename);

    AtomCollection getAtomCollection();

    SpinTypeCollection getSpinTypeCollection() const;
    ParticleCollection getParticleCollection(unsigned long k, unsigned long m) const;
    ElectronCollection getElectronCollection(unsigned long k, unsigned long m) const;
    ElectronCollections getElectronCollections(unsigned long k) const;

    unsigned long getNumberOfMaxima(unsigned long k, unsigned long m) const;
    double getNegativeLogarithmizedProbabilityDensity(unsigned long k, unsigned long m) const;

private:
    void countSubstructures();
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

class PathFileImporter : public Importer{

};

class WfFileImporter : public Importer{
public:
    WfFileImporter(const std::string& filename);

private:

};

#endif //AMOLQCGUI_IMPORTER_H
