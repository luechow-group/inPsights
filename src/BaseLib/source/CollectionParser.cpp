//
// Created by Michael Heuer on 08.12.17.
//

#include "CollectionParser.h"
#include <fstream>
#include <iostream>
#include "ElementInfo.h"

using json = nlohmann::json;

nlohmann::json CollectionParser::particleCollectionToJson(const ParticleCollection &particleCollection) {
    json j;
    auto particleCoordinatesArray = json::array();
    for (int i = 0; i < particleCollection.numberOfParticles(); ++i) {
        auto vec = particleCollection[i].position();
        particleCoordinatesArray.push_back(json::array({vec[0],vec[1],vec[2]}).dump());
    }

    j["coordinates"] = particleCoordinatesArray;
    return j;
}

nlohmann::json CollectionParser::atomCollectionToJson(const AtomCollection &atomCollection) {

    json j;
    auto elementTypeCollection = static_cast<ElementTypeCollection>(atomCollection);

    auto elementSymbolArray = json::array();
    auto particleCoordinatesArray = json::array();
    for (int i = 0; i < atomCollection.numberOfParticles(); ++i) {

        auto vec = atomCollection[i].position();
        auto elemetSymbol = Elements::ElementInfo::symbol(elementTypeCollection.elementType(i));

        elementSymbolArray.push_back(elemetSymbol);
        particleCoordinatesArray.push_back(json::array({vec[0],vec[1],vec[2]}).dump());
    }

    j["elements"] = elementSymbolArray;
    j["coordinates"] = particleCoordinatesArray;

    return j;
}

nlohmann::json CollectionParser::electronCollectionToJson(const ElectronCollection &electronCollection) {

    json j;
    auto spinTypeCollection = static_cast<SpinTypeCollection>(electronCollection);

    auto spinTypeArray = json::array();
    auto particleCoordinatesArray = json::array();
    for (int i = 0; i < electronCollection.numberOfParticles(); ++i) {

        auto vec = electronCollection[i].position();
        spinTypeArray.push_back(spinTypeCollection.spinType(i));
        particleCoordinatesArray.push_back(json::array({vec[0],vec[1],vec[2]}).dump());
    }

    j["spins"] = spinTypeArray.dump();
    j["coordinates"] = particleCoordinatesArray;

    return j;
}

void CollectionParser::writeJSON(const json& json, const std::string& filename) {
    std::ofstream ofstream(filename);
    ofstream << json.dump(4) << std::endl;
    ofstream.close();
}

//ParticleCollection CollectionParser::readParticleCollection(std::string &filename) {
//    return ParticleCollection();}

/*
bool is_empty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}
 */