//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include "ReferenceData.h"

int main(int argc, char *argv[]) {

    ReferenceSampleMapping mapping;

    RawDataReader reader(mapping);
    reader.read("raw1.bin");


    for (const auto& i : mapping.map){
        std::cout << i.first.maximum_ << std::endl;
    }
}