//
// Created by Michael Heuer on 27.08.18.
//

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
