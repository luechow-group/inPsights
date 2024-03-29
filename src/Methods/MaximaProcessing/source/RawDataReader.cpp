// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <RawDataReader.h>
#include "Maximum.h"
#include "ParticlesVector.h"
#include "ElectronSelection.h"
#include "MaximaProcessingSettings.h"
#include <iomanip>
#include <fstream>
#include <random>

RawDataReader::RawDataReader(
        Cluster& references,
        std::vector<Sample>& samples, int recordDelimiterLength)
        :
        BinaryFileReader(recordDelimiterLength),
        maxima_(references),
        samples_(samples),
        id_(0)
        {}


void RawDataReader::read(const std::string &basename) {
    read(basename, std::numeric_limits<size_t>::max(), false); // read all samples
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

void RawDataReader::readSamplesAndMaxima(std::ifstream &input, int fileLength, size_t numberOfSamples, bool shuffle) {
    auto ne = spins_.numberOfEntities();

    while (checkEOF(input, fileLength) && id_ < numberOfSamples) {
        auto serializedData = readVectorXd(input, size_t(ne)*7+1, 1);
        assert(serializedData.allFinite() && "The imported data must contain only finite values (non NaN or +/- Inf).");

        // create sample
        auto sample = serializedData.segment(0,3*ne);
        auto kineticEnergies = serializedData.segment(3*ne, ne);
        auto s = Sample(ElectronsVector(PositionsVector(sample), spins_),kineticEnergies);

        // create reference
        auto maximum =  serializedData.segment(4*ne, 3*ne);
        auto value = serializedData[7*ne];
        auto r = Maximum(atoms_, value, ElectronsVector(PositionsVector(maximum), spins_), id_);


        if(MaximaProcessing::settings.deleteCoreElectrons()) {
            auto [sValenceOnly, rValenceOnly] = removeNonValenceElectrons(s, r);
            samples_.emplace_back(std::move(sValenceOnly));
            maxima_.emplace_back(Cluster(std::move(rValenceOnly)));
        } else {
            samples_.emplace_back(std::move(s));
            maxima_.emplace_back(std::move(r));
        }

        id_++;
    }

    if (shuffle) {
        // Random shuffled vector to remove method artifacts
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm(ne);
        perm.setIdentity();
        for (auto &ref: maxima_) {
            std::shuffle(perm.indices().data(), perm.indices().data() + perm.indices().size(),
                         std::mt19937(MaximaProcessing::settings.seed()));
            ref.permuteAll(perm, samples_);
        }
    }
}

void RawDataReader::read(const std::string &basename, size_t numberOfSamples, bool shuffle){

    unsigned fileCounter = 0;
    int fileLength;

    auto filename = getFilename(basename, fileCounter);

    std::ifstream input(filename.c_str(), std::ios::binary);
    if( input.good() ) {
        fileLength = static_cast<int>(getFileLength(input));
    }
    else throw std::runtime_error("Could not open file '" + filename + "'");

    readHeader(input);

    auto samplesLimit = (numberOfSamples == 0)? std::numeric_limits<size_t>::max() : numberOfSamples;

    while(input.good() && maxima_.size() < samplesLimit) {
            readSamplesAndMaxima(input, fileLength, samplesLimit, shuffle);
            input.close();

            fileCounter++;
            filename = getFilename(basename, fileCounter);
            input = std::ifstream(filename.c_str(), std::ios::binary);

            fileLength = static_cast<int>(getFileLength(input));
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
    assert(nAlpha >= 0 && nElectrons > 0);
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

std::pair<Sample, Maximum>
RawDataReader::removeNonValenceElectrons(const Sample &sample, const Maximum &reference) {

    auto nonValenceIndices = ElectronSelection::getNonValenceIndices(reference.maximum(), atoms_);

    ElectronsVector newMaximum, newSample;
    Eigen::VectorXd newKineticEnergies(reference.maximum().numberOfEntities()-nonValenceIndices.size());

    long count = 0;
    for (long i = 0; i < reference.maximum().numberOfEntities(); ++i) {
        if(auto it = std::find(nonValenceIndices.begin(),nonValenceIndices.end(),i); it == nonValenceIndices.end()){
            newMaximum.append(reference.maximum()[i]);
            newSample.append(sample.sample_[i]);
            newKineticEnergies[count] = sample.kineticEnergies_[i];
            ++count;
        }
    }

    return {Sample (newSample, newKineticEnergies), Maximum(atoms_, reference.value(), newMaximum, reference.ownId())};
}
