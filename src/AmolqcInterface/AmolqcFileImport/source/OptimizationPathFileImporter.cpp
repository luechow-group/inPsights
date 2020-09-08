// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <OptimizationPathFileImporter.h>
#include <ParticlesVectorCollection.h>
#include <ElementInfo.h>

OptimizationPathFileImporter::OptimizationPathFileImporter(const std::string &filename)
        : AmolqcImporter(filename),atoms_() {

    numberOfNuclei_ = std::stoul(split(getLine(0))[0]);

    beginOfElectronPositionBlocks_ = numberOfNuclei_+2;

    numberOfAlphaElectrons_ = std::stoul(split(getLine(beginOfElectronPositionBlocks_-1))[5]);
    numberOfBetaElectrons_ = std::stoul(split(getLine(beginOfElectronPositionBlocks_-1))[7]);
    numberOfElectrons_= numberOfAlphaElectrons_+numberOfBetaElectrons_;


    for (unsigned i = 1; i <= numberOfNuclei_; ++i) {
        std::vector<std::string> lineElements = split(getLine(i));
        Element elementType = Elements::ElementInfo::elementTypeFromSymbol(lineElements[1]);
        double x = std::stod(lineElements[2]);
        double y = std::stod(lineElements[3]);
        double z = std::stod(lineElements[4]);

        atoms_.append({elementType,{x,y,z}});
    }

    unsigned long blockLength = numberOfElectrons_+2;
    substructuresData_ = countSubstructures(beginOfElectronPositionBlocks_,blockLength);
    numberOfPaths_ = substructuresData_.size();
}

unsigned long OptimizationPathFileImporter::calculateLine(unsigned long k, unsigned long m) const {
    assert( k > 0 && "k value must be greater than zero");
    assert( m > 0 && "m value must be greater than zero");
    assert( k <= numberOfPaths_ && "k value must be smaller than kmax");
    assert( m <= substructuresData_[k-1].numberOfSubstructures_ && "m value must be smaller than mmax");

    unsigned long start = substructuresData_[k-1].startingLine_;
    unsigned long linesToSkip = (m-1)*(numberOfElectrons_+2);
    return start + linesToSkip;
}

PositionsVector OptimizationPathFileImporter::getPositionsVector(unsigned long k, unsigned long m) const {
    unsigned long startLine = calculateLine(k,m)+2;
    return AmolqcImporter::importPositionsVectorBlock(startLine, 0, numberOfElectrons_);
}

ElectronsVectorCollection OptimizationPathFileImporter::getPath(unsigned long k) const {
    unsigned long numberOfSubstructures = substructuresData_[k-1].numberOfSubstructures_;
    PositionsVectorCollection positionsVectorCollection;
    for (unsigned long m = 1; m <= numberOfSubstructures; ++m) {
        positionsVectorCollection.append(this->getPositionsVector(k, m));
    }
    return ElectronsVectorCollection(positionsVectorCollection,
                               SpinTypesVector(numberOfAlphaElectrons_, numberOfBetaElectrons_));
}

unsigned long OptimizationPathFileImporter::getNumberOfPaths() const {
    return numberOfPaths_;
}

AtomsVector OptimizationPathFileImporter::getAtomsVector() const {
    return atoms_;
}
