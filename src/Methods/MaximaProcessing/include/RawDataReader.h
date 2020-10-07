// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_RAWDATAREADER_H
#define INPSIGHTS_RAWDATAREADER_H

#include <BinaryFileReader.h>
#include "Maximum.h"
#include "Sample.h"

class RawDataReader : public BinaryFileReader{
public:
    explicit RawDataReader(
            Cluster& maxima,
            std::vector<Sample>& samples,
            int recordDelimiterLength = 4);

    void read(const std::string& basename) override;
    
    void read(const std::string& basename, size_t numberOfSamples);

    std::pair<Sample, Maximum> removeNonValenceElectrons(const Sample &sample, const Maximum &reference);

    AtomsVector getAtoms() const;
    
private:
    void readHeader(std::ifstream& input);
    void readAtomsHeader(std::ifstream &input);
    void readElectronsHeader(std::ifstream &input);
    void readSamplesAndMaxima(std::ifstream &input, int fileLength, size_t numberOfSamples);
    std::string zeroPadNumber(int num);
    long long int getFileLength(std::ifstream &input) const;
    std::string getFilename(const std::string &basename, unsigned fileCounter);

    AtomsVector atoms_;
    SpinTypesVector spins_;
    Cluster& maxima_;
    std::vector<Sample>& samples_;
    size_t id_;
};

#endif //INPSIGHTS_RAWDATAREADER_H
