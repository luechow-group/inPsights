//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include "ReferenceSampleMapping.h"

int main(int argc, char *argv[]) {

    ReferenceSampleMapping mapping;

    RawDataReader reader(mapping);
    reader.read("raw1.bin");


    Reference r;
    for (auto i : mapping.map){
        r=Reference(i.first.maximum_,i.first.negLogSqrdProbabilityDensity_);
    }

    auto res = mapping.map.find(r);
    if(res == mapping.map.end())
        std::cout << "END" << std::endl;
    else{
        std::cout << (*res).first.maximum_.typesVector() << std::endl;

    }




    //auto b = (*mapping.map.begin()).first < (*(mapping.map.begin())).first;
    //std::cout << b << std::endl;// <mapping.map

}