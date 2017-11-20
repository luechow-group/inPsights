//
// Created by Michael Heuer on 26.01.17.
//

#include <iostream>

#include "ElectronicWaveFunction.h"
#include "WfFileImporter.h"

int main(int argc, char const *argv[]) {
    std::cout << "Hello world" << std::endl;


    WfFileImporter importer("Ethane-em-5.wf");
    std::cout << importer.getBasisSet() << std::endl;

}
