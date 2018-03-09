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
    j["type"] = "ParticlesVector";
    j["coordinates"] = particleCoordinatesArray;
    return j;
}

nlohmann::json CollectionParser::atomsVectorToJson(const AtomsVector &atomsVector) {
    nlohmann::json j;
    const auto &elementTypesVector = atomsVector.elementTypesVector();

    auto elementSymbolArray = nlohmann::json::array();
    auto particleCoordinatesArray = nlohmann::json::array();
    for (int i = 0; i < atomsVector.numberOfEntities(); ++i) {

        auto vec = atomsVector[i].position();
        auto elemetSymbol = Elements::ElementInfo::symbol(elementTypesVector[i]);

        elementSymbolArray.push_back(elemetSymbol);
        particleCoordinatesArray.push_back(nlohmann::json::array({vec[0],vec[1],vec[2]}));
    }

    j["type"] = "AtomsVector";
    j["elements"] = elementSymbolArray;
    j["coordinates"] = particleCoordinatesArray;

    return j;
}

nlohmann::json CollectionParser::electronsVectorToJson(const ElectronsVector &electronsVector) {
    nlohmann::json j;
    const auto &spinTypesVector = electronsVector.spinTypesVector();

    auto spinTypeArray = nlohmann::json::array();
    auto particleCoordinatesArray = nlohmann::json::array();
    for (int i = 0; i < electronsVector.numberOfEntities(); ++i) {

        auto vec = electronsVector[i].position();
        spinTypeArray.emplace_back((int)spinTypesVector[i]);
        particleCoordinatesArray.emplace_back(nlohmann::json::array({vec[0],vec[1],vec[2]}));
    }

    j["type"] = "ElectronsVector";
    j["spins"] = spinTypeArray;
    j["coordinates"] = particleCoordinatesArray;

    return j;
}


PositionsVector CollectionParser::positionsVectorFromJson(const std::string &filename) {
    auto j = readJSON(filename);
    assert(j["type"]== "ParticlesVector" && "File must be a ParticlesVector.");

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


AtomsVector CollectionParser::atomsVectorFromJson(const std::string &filename) {
    auto j = readJSON(filename);
    assert(j["type"]== "AtomsVector" && "File must be a AtomsVector.");
    assert(j["elements"].is_array());

    auto elementSymbols = j["elements"].get<std::vector<std::string>>();
    ElementTypesVector elementTypesVector;
    for(auto it = elementSymbols.begin(); it != elementSymbols.end(); ++it) {
        elementTypesVector.append(Elements::ElementInfo::elementTypeForSymbol(*it));
    }
    PositionsVector positionsVector = array2DToPositionsVector(j["coordinates"]);

    return AtomsVector(positionsVector, elementTypesVector);
}

ElectronsVector CollectionParser::electronsVectorFromJson(const std::string &filename) {
    auto j = readJSON(filename);
    assert(j["type"]== "ElectronsVector" && "File must be a ElectronsVector.");
    assert(j["spins"].is_array());

    auto spins = j["spins"].get<std::vector<int>>();
    SpinTypesVector spinTypesVector;
    for(auto it = spins.begin(); it != spins.end(); ++it) {
        spinTypesVector.append( (Spin::SpinType)(*it) );
    }
    PositionsVector positionsVector = array2DToPositionsVector(j["coordinates"]);

    return ElectronsVector(positionsVector, spinTypesVector);
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
