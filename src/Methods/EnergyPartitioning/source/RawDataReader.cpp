//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include "Reference.h"
#include "ParticlesVector.h"
#include <Logger.h>
#include <spdlog/spdlog.h>
#include <experimental/filesystem>
#include <iomanip>

namespace fs = std::experimental::filesystem;

RawDataReader::RawDataReader(
        std::vector<Reference>& references,
        std::vector<Sample>& samples, int recordDelimiterLength)
        :
        BinaryFileReader(recordDelimiterLength),
        references_(references),
        samples_(samples)
        {}


void RawDataReader::read(const std::string &basename) {
    read(basename, std::numeric_limits<size_t>::max()); // read all samples
}

void RawDataReader::readAtomsHeader(std::ifstream &input) {

    int nAtoms = readInt(input);
    assert(nAtoms >=0);
    auto numberOfAtoms = static_cast<unsigned>(nAtoms);

    Eigen::VectorXi atomicNumbers(numberOfAtoms);
    Eigen::VectorXd atomCoords(numberOfAtoms*3);

    atomicNumbers = readVectorXi(input,size_t(nAtoms),1);
    atomCoords = readVectorXd(input, size_t(nAtoms),3);

    std::vector<Element> elementTypes(numberOfAtoms);
    for (size_t i = 0; i < numberOfAtoms; ++i) {
        elementTypes[i] = Elements::elementFromInt(atomicNumbers[i]);
    }

    atoms_ = AtomsVector(PositionsVector(atomCoords), elementTypes);
}

void RawDataReader::readSamplesAndMaxima(std::ifstream &input, int fileLength, size_t numberOfSamples) {
    auto ne = spins_.numberOfEntities();

    size_t id = 0;

    while (checkEOF(input, fileLength) && id < numberOfSamples) {

        auto serializedData = readVectorXd(input, size_t(ne)*7+1, 1);

        // create sample
        auto sample = serializedData.segment(0,3*ne);
        auto kineticEnergies = serializedData.segment(3*ne, ne);
        auto s = Sample(ElectronsVector(PositionsVector(sample), spins_),kineticEnergies);
        samples_.emplace_back(s);

        // create reference
        auto maximum =  serializedData.segment(4*ne, 3*ne);
        auto value = serializedData[7*ne];
        auto r = Reference(value, ElectronsVector(PositionsVector(maximum), spins_), id);
        references_.emplace_back(r);

        id++;
    }
}

void RawDataReader::read(const std::string &basename, size_t numberOfSamples){

    std::ifstream input;
    unsigned fileCounter = 0;
    int fileLength;

    auto filename = getFilename(basename, fileCounter);
    if( fs::exists(filename) ) {
        input = std::ifstream(filename.c_str(), std::ios::binary);
        fileLength = static_cast<int>(getFileLength(input));
    }
    else throw std::runtime_error("Could not open file '" + basename + "'");

    if( input.good() ) readHeader(input);
    else throw std::runtime_error("Could not open file '" + basename + "'");

    while(fs::exists(filename) ) {
        if( input.good() ) {
            readSamplesAndMaxima(input, fileLength, numberOfSamples);

            fileCounter++;
            filename = getFilename(basename,fileCounter);
            input = std::ifstream(filename.c_str(), std::ios::binary);
            fileLength = static_cast<int>(getFileLength(input));
        }
        else throw std::runtime_error("Could not open file '" + basename + "'");
    }
}

std::string RawDataReader::getFilename(const std::string &basename, unsigned fileCounter) {
    return basename + "-" + zeroPadNumber(fileCounter) + ".bin";
}

long long int RawDataReader::getFileLength(std::ifstream &input) const {// get length of file:
    input.seekg(0, std::ios_base::end);
    long long int totalLength = input.tellg();
    input.seekg(0, std::ios_base::beg);
    return totalLength;
}

void RawDataReader::readHeader(std::ifstream &input){
    readAtomsHeader(input);
    readElectronsHeader(input);
}

void RawDataReader::readElectronsHeader(std::ifstream &input) {// electrons header
    int nElectrons = readInt(input);
    int nAlpha = readInt(input);
    assert(nAlpha >=0 && nElectrons > 0);
    spins_ = SpinTypesVector(static_cast<unsigned>(nAlpha), static_cast<unsigned>(nElectrons - nAlpha));
}

AtomsVector RawDataReader::getAtoms() const {
    return atoms_;
}

std::string RawDataReader::zeroPadNumber(int num) {
    std::ostringstream ss;
    ss << std::setw(2) << std::setfill('0') << num;
    return ss.str();
}
