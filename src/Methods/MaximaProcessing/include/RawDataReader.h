/* Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef INPSIGHTS_RAWDATAREADER_H
#define INPSIGHTS_RAWDATAREADER_H

#include <BinaryFileReader.h>
#include "Reference.h"
#include "Sample.h"

class RawDataReader : public BinaryFileReader{
public:
    explicit RawDataReader(
            Group& maxima,
            std::vector<Sample>& samples,
            int recordDelimiterLength = 4);

    void read(const std::string& basename) override;
    
    void read(const std::string& basename, size_t numberOfSamples);

    std::pair<Sample, Reference> removeNonValenceElectrons(const Sample &sample, const Reference &reference);

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
    Group& maxima_;
    std::vector<Sample>& samples_;
    size_t id_;
};

#endif //INPSIGHTS_RAWDATAREADER_H
