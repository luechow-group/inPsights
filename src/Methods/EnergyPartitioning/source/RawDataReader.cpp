//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>


RawDataReader::RawDataReader(int recordDelimiterLength)
        : BinaryFileReader(recordDelimiterLength)
{}

bool RawDataReader::read(const std::string &fileName){
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

