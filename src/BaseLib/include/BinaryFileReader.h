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

#ifndef INPSIGHTS_BINARYFILEREADER_H
#define INPSIGHTS_BINARYFILEREADER_H

#include <ParticlesVector.h>
#include <fstream>

class BinaryFileReader {
public:

    explicit BinaryFileReader(int recordDelimiterLength = 4);

    virtual void read(const std::string &fileName) = 0;

    bool checkEOF(std::ifstream &input, long long int totalLength) const;

    int readInt(std::ifstream &input) const;

    double readDouble(std::ifstream &input) const;

    Eigen::VectorXd readVectorXd(std::ifstream &input, size_t numberOfEntities, size_t entityLength = 1) const;

    Eigen::VectorXi readVectorXi(std::ifstream &input, size_t numberOfEntities, size_t entityLength = 1) const;

private:
    int recordDelimiterLength_ = 4;
};

#endif //INPSIGHTS_BINARYFILEREADER_H
