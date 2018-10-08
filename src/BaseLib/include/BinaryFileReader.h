//
// Created by Michael Heuer on 27.08.18.
//

#ifndef AMOLQCPP_BINARYFILEREADER_H
#define AMOLQCPP_BINARYFILEREADER_H

#include <ParticlesVector.h>
#include <vector>
#include <iostream>
#include <fstream>


class BinaryFileReader{
public:

    explicit BinaryFileReader(int recordDelimiterLength = 4)
    : recordDelimiterLength_(recordDelimiterLength){};


    virtual void read(const std::string& fileName) = 0;

    bool checkEOF(std::ifstream &input, long long int totalLength){
        return (totalLength-input.tellg()) >= recordDelimiterLength_;
    }


    int readInt(std::ifstream &input) const {
        int value;
        input.seekg(recordDelimiterLength_, std::ios::cur);
        input.read((char*) &value, sizeof(int));
        input.seekg(recordDelimiterLength_, std::ios::cur);
        return value;
    }

    double readDouble(std::ifstream &input) const {
        double value;
        input.seekg(recordDelimiterLength_, std::ios_base::cur);
        input.read((char*) &value, sizeof(double) );
        input.seekg(recordDelimiterLength_, std::ios_base::cur);
        return value;
    }

    Eigen::VectorXd readVectorXd(std::ifstream &input,  size_t numberOfEntities, size_t entityLength = 1) const {
        double coords[entityLength*numberOfEntities];
        input.seekg(recordDelimiterLength_, std::ios_base::cur);
        input.read((char*) &coords, entityLength*numberOfEntities*sizeof(double) );
        input.seekg(recordDelimiterLength_, std::ios_base::cur);
        return Eigen::Map<Eigen::VectorXd>(coords, entityLength*numberOfEntities);
    }

private:
    int recordDelimiterLength_ = 4;
};

#endif //AMOLQCPP_BINARYFILEREADER_H
