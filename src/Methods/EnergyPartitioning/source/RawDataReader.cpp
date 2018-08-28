//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include "ReferenceData.h"
#include "ParticlesVector.h"

RawDataReader::RawDataReader(ReferenceSampleMapping& mapping, int recordDelimiterLength)
        : BinaryFileReader(recordDelimiterLength), mapping_(mapping) {}

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

        unsigned numberOfAlphaElectrons, numberOfBetaElectrons;;

        assert(nAlpha >=0 && nElectrons > 0);
        numberOfAlphaElectrons = static_cast<unsigned>(nAlpha);
        numberOfBetaElectrons = static_cast<unsigned>(nElectrons-nAlpha);
        auto spins = SpinTypesVector(numberOfAlphaElectrons,numberOfBetaElectrons);


        std::cout << nElectrons << " " << nAlpha << std::endl; // N alpha electrons, core positions

        while (checkEOF(input,totalLength)) {



            auto s = Sample(
                    ElectronsVector(PositionsVector(readVectorXd(input, size_t(nElectrons),3)), spins),
                    readVectorXd(input, size_t(nElectrons)));

            auto r = Reference(
                    ElectronsVector(PositionsVector(readVectorXd(input, size_t(nElectrons),3)), spins),
                    readDouble(input));
            //std::cout << readVectorXd(input, size_t(nElectrons), 3).transpose() << std::endl << std::endl;
            //std::cout << readVectorXd(input, size_t(nElectrons)).transpose() << std::endl << std::endl;
            //std::cout << readVectorXd(input, size_t(nElectrons), 3).transpose() << std::endl << std::endl;
            //std::cout << readDouble(input) << std::endl;

            mapping_.map.insert(RefSamplePair(r,s));

        }
        return true;
    }
    else
    {
        std::cout << "input not good"<< std::endl;
        return false;
    }


}

