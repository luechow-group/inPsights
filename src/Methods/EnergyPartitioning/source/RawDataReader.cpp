//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include "Reference.h"
#include "ParticlesVector.h"
#include <Logger.h>
#include <spdlog/spdlog.h>

RawDataReader::RawDataReader(
        std::vector<Reference>& references,
        std::vector<Sample>& samples, int recordDelimiterLength)
        :
        BinaryFileReader(recordDelimiterLength),
        references_(references),
        samples_(samples)
        {}

void RawDataReader::readAtoms(std::ifstream &input) {
    int nAtoms = readInt(input);
    assert(nAtoms >=0);
    auto numberOfAtoms = static_cast<unsigned>(nAtoms);

    Eigen::VectorXi atomicNumbers(numberOfAtoms);
    Eigen::VectorXd atomCoords(numberOfAtoms*3);

    atomicNumbers = readVectorXi(input,size_t(nAtoms),1);
    atomCoords = readVectorXd(input, size_t(nAtoms),3);

    std::vector<Element> elementTypes(numberOfAtoms);
    for (int i = 0; i < numberOfAtoms; ++i) {
        elementTypes[i] = Elements::elementFromInt(atomicNumbers[i]);
    }

    atoms_ = AtomsVector(PositionsVector(atomCoords), elementTypes);
}

void RawDataReader::readElectrons(std::ifstream &input, int totalLength) {
    int nElectrons = readInt(input);
    int nAlpha = readInt(input);

    unsigned numberOfAlphaElectrons, numberOfBetaElectrons;

    assert(nAlpha >=0 && nElectrons > 0);
    numberOfAlphaElectrons = static_cast<unsigned>(nAlpha);
    numberOfBetaElectrons = static_cast<unsigned>(nElectrons-nAlpha);
    auto spins = SpinTypesVector(numberOfAlphaElectrons,numberOfBetaElectrons);

    size_t id = 0;

    while (checkEOF(input, totalLength)) {

        // don't move read methods into constructor as this messes up the ifstream stride
        auto sample = readVectorXd(input, size_t(nElectrons),3);
        auto kineticEnergies = readVectorXd(input, size_t(nElectrons));
        auto s = Sample(ElectronsVector(PositionsVector(sample), spins),kineticEnergies);
        samples_.emplace_back(s);

        auto maximum = readVectorXd(input, size_t(nElectrons),3);
        auto value = readDouble(input);
        auto r = Reference(value, ElectronsVector(PositionsVector(maximum), spins), id);
        references_.emplace_back(r);

        id++;
    }
}

void RawDataReader::read(const std::string &fileName){
    // open file in binary mode
    std::ifstream input(fileName.c_str(), std::ios::binary);

    if( input.good() ) {
        // get length of file:
        input.seekg (0, std::ifstream::end);
        long long int totalLength = input.tellg();
        input.seekg (0, std::ifstream::beg);
        spdlog::get(Logger::name)->info("file length {}", totalLength);

        readAtoms(input);

        readElectrons(input, totalLength);
    } else
        throw std::runtime_error("Could not open file '" + fileName + "'");


}

AtomsVector RawDataReader::getAtoms() const {
    return atoms_;
}
