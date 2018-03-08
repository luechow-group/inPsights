//
// Created by Michael Heuer on 08.12.17.
//

#include "CollectionParser.h"
#include <fstream>
#include "ElementInfo.h"

nlohmann::json CollectionParser::positionCollectionToJson(const PositionCollection &positionCollection) {
    nlohmann::json j;
    auto particleCoordinatesArray = nlohmann::json::array();
    for (int i = 0; i < positionCollection.numberOfEntities(); ++i) {
        auto vec = positionCollection[i];
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
    const auto &spinTypeCollection = electronCollection.spinTypeCollection();

    auto spinTypeArray = nlohmann::json::array();
    auto particleCoordinatesArray = nlohmann::json::array();
    for (int i = 0; i < electronCollection.numberOfEntities(); ++i) {

        auto vec = electronCollection[i].position();
        spinTypeArray.emplace_back((int)spinTypeCollection[i]);
        particleCoordinatesArray.emplace_back(nlohmann::json::array({vec[0],vec[1],vec[2]}));
    }

    j["type"] = "ElectronCollection";
    j["spins"] = spinTypeArray;
    j["coordinates"] = particleCoordinatesArray;

    return j;
}


PositionCollection CollectionParser::positionCollectionFromJson(const std::string &filename) {
    auto j = readJSON(filename);
    assert(j["type"]== "ParticleCollection" && "File must be a ParticleCollection.");

    return array2DToPositionCollection(j["coordinates"]);
}

PositionCollection CollectionParser::array2DToPositionCollection(const nlohmann::json &coordinates) {
    auto particles = coordinates.get<std::vector<std::vector<double>>>();
    PositionCollection positionCollection;
    for(auto it = particles.begin(); it != particles.end(); ++it) {
        assert((*it).size() == 3);
        positionCollection.append(Eigen::Vector3d((*it)[0],(*it)[1],(*it)[2]));
    }
    return positionCollection;
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
    PositionCollection positionCollection = array2DToPositionCollection(j["coordinates"]);

    return AtomCollection(positionCollection, elementTypeCollection);
}

ElectronCollection CollectionParser::electronCollectionFromJson(const std::string &filename) {
    auto j = readJSON(filename);
    assert(j["type"]== "ElectronCollection" && "File must be a ElectronCollection.");
    assert(j["spins"].is_array());

    auto spins = j["spins"].get<std::vector<int>>();
    SpinTypeCollection spinTypeCollection;
    for(auto it = spins.begin(); it != spins.end(); ++it) {
        spinTypeCollection.append( (Spin::SpinType)(*it) );
    }
    PositionCollection positionCollection = array2DToPositionCollection(j["coordinates"]);

    return ElectronCollection(positionCollection, spinTypeCollection);
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
