//
// Created by Michael Heuer on 06.11.17.
//

#include "AmolqcFileImport/RefFileImporter.h"
#include "ElementInfo.h"

RefFileImporter::RefFileImporter(const std::string &filename)
        : AmolqcImporter(filename) {

    numberOfNuclei_ = std::stoul(split(getLine(0))[0]);
    numberOfElectrons_ = std::stoul(split(getLine(numberOfNuclei_+2))[10]);
    numberOfAlphaElectrons_ = std::stoul(split(getLine(numberOfNuclei_+2))[8]);
    numberOfBetaElectrons_ = numberOfElectrons_-numberOfAlphaElectrons_;
    std::vector<std::string> summaryLineElements = split(getLine(numberOfNuclei_+1));

    numberOfSuperstructures_ = std::stoul(summaryLineElements[0]);
    totalNumberOfMaxima_ = std::stoul(summaryLineElements[2]);
    maximalNumberOfSubstructures = std::stoul(summaryLineElements[3]);


    unsigned long startLineIdx = numberOfNuclei_+2;
    unsigned long blockLength = numberOfElectrons_+2;
    substructuresData_ = AmolqcImporter::countSubstructures(startLineIdx,blockLength);
}

AtomsVector RefFileImporter::getAtomsVector() {
    AtomsVector atomsVector;

    for (unsigned i = 1; i <= numberOfNuclei_; ++i) {
        std::vector<std::string> lineElements = split(getLine(i));
        Element elementType = Elements::ElementInfo::elementTypeFromSymbol(lineElements[1]);
        double x = std::stod(lineElements[2]);
        double y = std::stod(lineElements[3]);
        double z = std::stod(lineElements[4]);
        atomsVector.append({elementType,{x,y,z}});
    }
    return atomsVector;
}

unsigned long RefFileImporter::calculateLine(unsigned long k, unsigned long m) const {
    assert( k > 0 && "k value must be greater than zero");
    assert( m > 0 && "m value must be greater than zero");
    assert( k <= numberOfSuperstructures_ && "k value must be smaller than kmax");
    assert( m <= substructuresData_[k-1].numberOfSubstructures_  && "m value must be smaller than or equal to mmax");

    unsigned long start = substructuresData_[k-1].startingLine_;
    unsigned long linesToSkip = (m-1)*(numberOfElectrons_+2);
    return start + linesToSkip;
}

PositionsVector RefFileImporter::getPositionsVector(unsigned long k, unsigned long m) const {
    unsigned long startLine = calculateLine(k,m)+numberOfLinesAboveCoordinatesBlock;
    return importPositionsVectorBlock(startLine, 0, numberOfElectrons_);
}

unsigned long RefFileImporter::getNumberOfMaxima(unsigned long k, unsigned long m) const {
    unsigned long startLine = calculateLine(k,m)+numberOfLinesAboveCoordinatesBlock;
    std::vector<std::string> lineElements = split(getLine(startLine));
    return std::stoul(lineElements[6]);
}

double RefFileImporter::getNegativeLogarithmizedProbabilityDensity(unsigned long k, unsigned long m) const {
    unsigned long startLine = calculateLine(k,m)+numberOfLinesAboveCoordinatesBlock;
    std::vector<std::string> lineElements = split(getLine(startLine));
    return std::stod(lineElements[4]);
}

SpinTypesVector RefFileImporter::getSpinTypesVector() const {
    return AmolqcImporter::getSpinTypesVector(numberOfAlphaElectrons_, numberOfBetaElectrons_);
}

ElectronsVector RefFileImporter::getMaximaStructure(unsigned long k, unsigned long m) const {
    return ElectronsVector(this->getPositionsVector(k,m), this->getSpinTypesVector());
}

ElectronsVectorCollection RefFileImporter::getAllSubstructures(unsigned long k) const {
    assert( k > 0 && "k value must be greater than zero");
    unsigned long numberOfSubstructures = substructuresData_[k-1].numberOfSubstructures_;
    PositionsVectorCollection positionsVectorCollection;
    for (unsigned long m = 1; m <= numberOfSubstructures; ++m) {
        positionsVectorCollection.append(this->getPositionsVector(k,m));
    }
    return ElectronsVectorCollection(positionsVectorCollection,this->getSpinTypesVector());
}

unsigned long RefFileImporter::numberOfSuperstructures() {
    return numberOfSuperstructures_;
}
