//
// Created by Michael Heuer on 08.12.17.
//

#include "CollectionParser.h"
#include <fstream>
#include "ElementInfo.h"

nlohmann::json CollectionParser::particleCollectionToJson(const ParticleCollection &particleCollection) {
    nlohmann::json j;
    auto particleCoordinatesArray = nlohmann::json::array();
    for (int i = 0; i < particleCollection.numberOfParticles(); ++i) {
        auto vec = particleCollection[i].position();
        particleCoordinatesArray.push_back(nlohmann::json::array({vec[0],vec[1],vec[2]}).dump());
    }
    j["type"] = "ParticleCollection";
    j["coordinates"] = particleCoordinatesArray;
    return j;
}

nlohmann::json CollectionParser::atomCollectionToJson(const AtomCollection &atomCollection) {
    nlohmann::json j;
    auto elementTypeCollection = static_cast<ElementTypeCollection>(atomCollection);

    auto elementSymbolArray = nlohmann::json::array();
    auto particleCoordinatesArray = nlohmann::json::array();
    for (int i = 0; i < atomCollection.numberOfParticles(); ++i) {

        auto vec = atomCollection[i].position();
        auto elemetSymbol = Elements::ElementInfo::symbol(elementTypeCollection.elementType(i));

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
    auto spinTypes = electronCollection.getSpinTypes();

    auto spinTypeArray = nlohmann::json::array();
    auto particleCoordinatesArray = nlohmann::json::array();
    for (int i = 0; i < electronCollection.numberOfParticles(); ++i) {

        auto vec = electronCollection[i].position();
        spinTypeArray.push_back(spinTypes(i));
        particleCoordinatesArray.push_back(nlohmann::json::array({vec[0],vec[1],vec[2]}));
    }

    j["type"] = "ElectronCollection";
    j["spins"] = spinTypeArray;
    j["coordinates"] = particleCoordinatesArray;

    return j;
}


ParticleCollection CollectionParser::particleCollectionFromJson(const std::string &filename) {
    auto j = readJSON(filename);
    assert(j["type"]== "ParticleCollection" && "File must be a ParticleCollection.");

    return array2DToParticleCollection(j["coordinates"]);
}

ParticleCollection CollectionParser::array2DToParticleCollection(nlohmann::json coordinates) {
    auto particles = coordinates.get<std::vector<std::vector<double>>>();
    ParticleCollection particleCollection;
    for(auto it = particles.begin(); it != particles.end(); ++it) {
        assert((*it).size() == 3);
        particleCollection.append(Particle((*it)[0],(*it)[1],(*it)[2]));
    }
    return particleCollection;
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
    ParticleCollection particleCollection = array2DToParticleCollection(j["coordinates"]);

    return AtomCollection(particleCollection, elementTypeCollection);
}

ElectronCollection CollectionParser::electronCollectionFromJson(const std::string &filename) {
    auto j = readJSON(filename);
    assert(j["type"]== "ElectronCollection" && "File must be a ElectronCollection.");
    assert(j["spins"].is_array());

    auto spins = j["spins"].get<std::vector<int>>();
    Eigen::VectorXi spinTypes(spins.size());
    //for(auto it = spins.begin(); it != spins.end(); ++it) {
    //    long i = std::distance(spins.begin(), it);
    for (int i = 0; i < spins.size(); ++i) {
        spinTypes[i] = spins.at(i);
    }
    ParticleCollection particleCollection = array2DToParticleCollection(j["coordinates"]);

    return ElectronCollection(particleCollection, spinTypes);
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
