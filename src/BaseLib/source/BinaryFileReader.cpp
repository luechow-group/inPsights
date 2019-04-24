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

#include <BinaryFileReader.h>

BinaryFileReader::BinaryFileReader(int recordDelimiterLength)
        : recordDelimiterLength_(recordDelimiterLength) {};

bool BinaryFileReader::checkEOF(std::ifstream &input, long long int totalLength) const {
    return (totalLength - input.tellg()) >= recordDelimiterLength_;
}

int BinaryFileReader::readInt(std::ifstream &input) const {
    int value;
    input.seekg(recordDelimiterLength_, std::ios::cur);
    input.read((char *) &value, sizeof(int));
    input.seekg(recordDelimiterLength_, std::ios::cur);
    return value;
}

double BinaryFileReader::readDouble(std::ifstream &input) const {
    double value;
    input.seekg(recordDelimiterLength_, std::ios_base::cur);
    input.read((char *) &value, sizeof(double));
    input.seekg(recordDelimiterLength_, std::ios_base::cur);
    return value;
}

Eigen::VectorXd BinaryFileReader::readVectorXd(
        std::ifstream &input,
        size_t numberOfEntities,
        size_t entityLength) const {

    double coords[entityLength * numberOfEntities];
    input.seekg(recordDelimiterLength_, std::ios_base::cur);
    input.read((char *) &coords, entityLength * numberOfEntities * sizeof(double));
    input.seekg(recordDelimiterLength_, std::ios_base::cur);
    return Eigen::Map<Eigen::VectorXd>(coords, entityLength * numberOfEntities);
}

Eigen::VectorXi BinaryFileReader::readVectorXi(
        std::ifstream &input,
        size_t numberOfEntities,
        size_t entityLength) const {

    int coords[entityLength * numberOfEntities];
    input.seekg(recordDelimiterLength_, std::ios_base::cur);
    input.read((char *) &coords, entityLength * numberOfEntities * sizeof(int));
    input.seekg(recordDelimiterLength_, std::ios_base::cur);
    return Eigen::Map<Eigen::VectorXi>(coords, entityLength * numberOfEntities);
}