//
// Created by Michael Heuer on 08.12.17.
//

#include "CollectionParser.h"
#include <fstream>
#include "ElementInfo.h"

nlohmann::json CollectionParser::positionsVectorToJson(const PositionsVector &positionsVector) {
    nlohmann::json j;
    auto particleCoordinatesArray = nlohmann::json::array();
    for (int i = 0; i < positionsVector.numberOfEntities(); ++i) {
        auto vec = positionsVector[i];
        particleCoordinatesArray.push_back(nlohmann::json::array({vec[0],vec[1],vec[2]}).dump());
    }
    j["type"] = "ParticleCollection";
    j["coordinates"] = particleCoordinatesArray;
    return j;
}

nlohmann::json CollectionParser::atomCollectionToJson(const AtomCollection &atomCollection) {
    nlohmann::json j;
    const auto &elementTypeCollection = atomCollection.elementTypeCollection();

    auto elementSymbolArray = nlohmann::json::array();
    auto particleCoordinatesArray = nlohmann::json::array();
    for (int i = 0; i < atomCollection.numberOfEntities(); ++i) {

        auto vec = atomCollection[i].position();
        auto elemetSymbol = Elements::ElementInfo::symbol(elementTypeCollection[i]);

        elementSymbolArray.push_back(elemetSymbol);
        particleCoordinatesArray.push_back(nlohmann::json::array({vec[0],vec[1],vec[2]}));
    }

    j["type"] = "AtomCollection";
    j["elements"] = elementSymbolArray;
    j["coordinates"] = particleCoordinatesArray;

    return j;
}

nlohmann::json CollectionParser::electronCollectionToJson(const ElectronCollection &electronCollection) {
    nlohmann::json j;
    const auto &spinTypesVector = electronCollection.spinTypesVector();

    auto spinTypeArray = nlohmann::json::array();
    auto particleCoordinatesArray = nlohmann::json::array();
    for (int i = 0; i < electronCollection.numberOfEntities(); ++i) {

        auto vec = electronCollection[i].position();
        spinTypeArray.emplace_back((int)spinTypesVector[i]);
        particleCoordinatesArray.emplace_back(nlohmann::json::array({vec[0],vec[1],vec[2]}));
    }

    j["type"] = "ElectronCollection";
    j["spins"] = spinTypeArray;
    j["coordinates"] = particleCoordinatesArray;

    return j;
}


PositionsVector CollectionParser::positionsVectorFromJson(const std::string &filename) {
    auto j = readJSON(filename);
    assert(j["type"]== "ParticleCollection" && "File must be a ParticleCollection.");

    return array2DToPositionsVector(j["coordinates"]);
}

PositionsVector CollectionParser::array2DToPositionsVector(const nlohmann::json &coordinates) {
    auto particles = coordinates.get<std::vector<std::vector<double>>>();
    PositionsVector positionsVector;
    for(auto it = particles.begin(); it != particles.end(); ++it) {
        assert((*it).size() == 3);
        positionsVector.append(Eigen::Vector3d((*it)[0],(*it)[1],(*it)[2]));
    }
    return positionsVector;
}


AtomCollection CollectionParser::atomCollectionFromJson(const std::string &filename) {
    auto j = readJSON(filename);
    assert(j["type"]== "AtomCollection" && "File must be a AtomCollection.");
    assert(j["elements"].is_array());

    auto elementSymbols = j["elements"].get<std::vector<std::string>>();
    ElementTypeCollection elementTypeCollection;
    for(auto it = elementSymbols.begin(); it != elementSymbols.end(); ++it) {
        elementTypeCollection.append(Elements::ElementInfo::elementTypeForSymbol(*it));
    }
    PositionsVector positionsVector = array2DToPositionsVector(j["coordinates"]);

    return AtomCollection(positionsVector, elementTypeCollection);
}

ElectronCollection CollectionParser::electronCollectionFromJson(const std::string &filename) {
    auto j = readJSON(filename);
    assert(j["type"]== "ElectronCollection" && "File must be a ElectronCollection.");
    assert(j["spins"].is_array());

    auto spins = j["spins"].get<std::vector<int>>();
    SpinTypesVector spinTypesVector;
    for(auto it = spins.begin(); it != spins.end(); ++it) {
        spinTypesVector.append( (Spin::SpinType)(*it) );
    }
    PositionsVector positionsVector = array2DToPositionsVector(j["coordinates"]);

    return ElectronCollection(positionsVector, spinTypesVector);
}

void CollectionParser::writeJSON(const nlohmann::json& json, const std::string& filename) {
    std::ofstream ofstream(filename);
    ofstream << json.dump(4) << std::endl;
    ofstream.close();
}

nlohmann::json CollectionParser::readJSON(const std::string &filename) {
    std::ifstream ifstream(filename);
    nlohmann::json j;
    j["filename"] = filename;
    ifstream >> j;

    return j;
}
