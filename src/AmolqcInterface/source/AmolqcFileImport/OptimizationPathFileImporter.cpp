//
// Created by Michael Heuer on 07.11.17.
//

#include "AmolqcFileImport/OptimizationPathFileImporter.h"
#include <ParticlesVectorCollection.h>

OptimizationPathFileImporter::OptimizationPathFileImporter(const std::string &filename,
                                                           unsigned long multiplicity)
        : AmolqcImporter(filename) {

    numberOfElectrons_= std::stoul(split(getLine(1))[0]);
    assert((numberOfElectrons_+(multiplicity-1))%2 == 0
           && "Even and oddness of the number of electrons must match to the specified multiplicity.");

    numberOfAlphaElectrons_ = (numberOfElectrons_+(multiplicity-1))/2;
    numberOfBetaElectrons_ = numberOfElectrons_-numberOfAlphaElectrons_;

    unsigned long blockLength = numberOfElectrons_+2;
    substructuresData_ = countSubstructures(0,blockLength);
    numberOfPaths_ = substructuresData_.size();
}

unsigned long OptimizationPathFileImporter::calculateLine(unsigned long k, unsigned long m) const {
    assert( k > 0 && "k value must be greater than zero");
    assert( m > 0 && "m value must be greater than zero");
    assert( k <= numberOfPaths_ && "k value must be smaller than kmax");
    assert( m <= substructuresData_[k].numberOfSubstructures_ && "m value must be smaller than mmax");

    unsigned long start = substructuresData_[k-1].startingLine_;
    unsigned long linesToSkip = (m-1)*(numberOfElectrons_+2);
    return start + linesToSkip;
}

PositionsVector OptimizationPathFileImporter::getPositionsVector(unsigned long k, unsigned long m) const {
    unsigned long startLine = calculateLine(k,m)+2;
    return AmolqcImporter::importPositionsVectorBlock(startLine, 0, numberOfElectrons_);
}

ElectronsVectorCollection OptimizationPathFileImporter::getPath(unsigned long k) const {
    unsigned long numberOfSubstructures = substructuresData_[k].numberOfSubstructures_;
    PositionsVectorCollection positionsVectorCollection;
    for (unsigned long m = 1; m <= numberOfSubstructures; ++m) {
        positionsVectorCollection.append(this->getPositionsVector(k, m));
    }
    return ElectronsVectorCollection(positionsVectorCollection,
                               this->getSpinTypesVector(numberOfAlphaElectrons_, numberOfBetaElectrons_));
}

unsigned long OptimizationPathFileImporter::getNumberOfPaths() const {
    return numberOfPaths_;
}
