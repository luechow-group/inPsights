//
// Created by heuer on 20.11.18.
//

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