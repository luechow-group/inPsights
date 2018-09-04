//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include "Reference.h"
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

        unsigned numberOfAlphaElectrons, numberOfBetaElectrons;

        assert(nAlpha >=0 && nElectrons > 0);
        numberOfAlphaElectrons = static_cast<unsigned>(nAlpha);
        numberOfBetaElectrons = static_cast<unsigned>(nElectrons-nAlpha);
        auto spins = SpinTypesVector(numberOfAlphaElectrons,numberOfBetaElectrons);

        while (checkEOF(input,totalLength)) {



            // don't move read methods into constructor as this messes up the ifstream stride
            auto sample = readVectorXd(input, size_t(nElectrons),3);
            auto kineticEnergies = readVectorXd(input, size_t(nElectrons));
            auto s = Sample(ElectronsVector(PositionsVector(sample), spins),kineticEnergies);

            auto maximum = readVectorXd(input, size_t(nElectrons),3);
            auto value = readDouble(input);
            auto r = Reference(ElectronsVector(PositionsVector(maximum), spins),value);

            mapping_.map.emplace(RefSamplePair(r,s));
            //https://stackoverflow.com/questions/21215214/how-to-sort-both-key-and-value-in-a-multimap

        }
        return true;
    }
    else
    {
        std::cout << "input not good"<< std::endl;
        return false;
    }



}

