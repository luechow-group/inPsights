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

    struct SamplePoint{
        ElectronsVector sample;
        Eigen::VectorXd kineticEnergies;
    };

    const int RECORD_DELIMITER_LENGTH = 4;

    bool read(const std::string& fileName)//, std::vector<RawDataElement>& points)
    {
        // clear the points
        //points.clear();

        // open file in binary mode
        std::ifstream input(fileName.c_str(), std::ios::binary);


        if( input.good() )
        {
            // get length of file:
            input.seekg (0, std::ifstream::end);
            long long int totalLength = input.tellg();
            input.seekg (0, std::ifstream::beg);
            std::cout << totalLength << std::endl; // in byte?


            int nElectrons = readInt(input);
            int nAlpha = readInt(input);
            std::cout << nElectrons << " " << nAlpha << std::endl; // N alpha electrons, core positions

            while (checkEOF(input,totalLength)) {
                std::cout << readVectorXd(input, size_t(nElectrons), 3).transpose() << std::endl << std::endl;
                std::cout << readVectorXd(input, size_t(nElectrons)).transpose() << std::endl << std::endl;
                std::cout << readVectorXd(input, size_t(nElectrons), 3).transpose() << std::endl << std::endl;
                std::cout << readDouble(input) << std::endl;

            }
            return true;
        }
        else
        {
            std::cout << "input not good"<< std::endl;
            return false;
        }
    }

    bool checkEOF(std::ifstream &input, long long int totalLength){
        return (totalLength-input.tellg()) >= RECORD_DELIMITER_LENGTH;
    }


    int readInt(std::ifstream &input) const {
        int value;
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        input.read((char*) &value, sizeof(int));
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios::cur);
        return value;
    }

    double readDouble(std::ifstream &input) const {
        double value;
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios_base::cur);
        input.read((char*) &value, sizeof(double) );
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios_base::cur);
        return value;
    }

    Eigen::VectorXd readVectorXd(std::ifstream &input,  size_t numberOfEntities, size_t entityLength = 1) const {
        double coords[entityLength*numberOfEntities];
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios_base::cur);
        input.read((char*) &coords, entityLength*numberOfEntities*sizeof(double) );
        input.seekg(RECORD_DELIMITER_LENGTH, std::ios_base::cur);
        return Eigen::Map<Eigen::VectorXd>(coords, entityLength*numberOfEntities);
    }

private:
};

#endif //AMOLQCPP_BINARYFILEREADER_H
