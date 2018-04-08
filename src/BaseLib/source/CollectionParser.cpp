//
// Created by Michael Heuer on 08.12.17.
//

#include "CollectionParser.h"
#include <fstream>
#include "ElementInfo.h"

nlohmann::json CollectionParser::positionsVectorToJson(const PositionsVector &positionsVector) {

    auto positions = nlohmann::json::array();
    for (int i = 0; i < positionsVector.numberOfEntities(); ++i) {
        positions.emplace_back(positionsVector[i][0]);
        positions.emplace_back(positionsVector[i][1]);
        positions.emplace_back(positionsVector[i][2]);
    }

    nlohmann::json j;
    j["PositionsVector"] = positions;
    return j;
}


nlohmann::json CollectionParser::elementTypesVectorToJson(const ElementTypesVector& elementTypesVector) {

    auto elementTypes = nlohmann::json::array();
    for (int i = 0; i < elementTypesVector.numberOfEntities(); ++i)
        elementTypes.emplace_back(Elements::ElementInfo::symbol(elementTypesVector[i]));

    nlohmann::json j;
    j["ElementTypesVector"] = elementTypes;
    return j;
}

nlohmann::json CollectionParser::spinTypesVectorToJson(const SpinTypesVector &spinTypesVector) {

    auto spinTypes = nlohmann::json::array();
    for (int i = 0; i < spinTypesVector.numberOfEntities(); i++){
        spinTypes.emplace_back((int)spinTypesVector[i]);
    }

    nlohmann::json j;
    j["SpinTypesVector"] = spinTypes;
    return j;
}


nlohmann::json CollectionParser::atomsVectorToJson(const AtomsVector &atomsVector) {
    nlohmann::json j;
    j["AtomsVector"]["ElementTypesVector"] = elementTypesVectorToJson(atomsVector.elementTypesVector())["ElementTypesVector"];
    j["AtomsVector"]["PositionsVector"] = positionsVectorToJson(atomsVector.positionsVector())["PositionsVector"];
    return j;
}

nlohmann::json CollectionParser::electronsVectorToJson(const ElectronsVector &electronsVector) {
    nlohmann::json j;
    j["ElectronsVector"]["SpinTypesVector"] = spinTypesVectorToJson(electronsVector.spinTypesVector())["SpinTypesVector"];
    j["ElectronsVector"]["PositionsVector"] = positionsVectorToJson(electronsVector.positionsVector())["PositionsVector"];
    return j;
}

nlohmann::json
CollectionParser::positionsVectorCollectionToJson(const PositionsVectorCollection &positionsVectorCollection) {

    auto positionsVectors = nlohmann::json::array();
    for (int i = 0; i < positionsVectorCollection.numberOfEntities(); i++) {
        positionsVectors.emplace_back(positionsVectorToJson(positionsVectorCollection[i])["PositionsVector"]);
    }

    nlohmann::json j;
    j["PositionsVectorCollection"] = positionsVectors;
    return j;
}

nlohmann::json CollectionParser::electronsVectorCollectionToJson(const ElectronsVectorCollection &electronsVectorCollection) {

    nlohmann::json j;
    j["ElectronsVectorCollection"]["SpinTypesVector"] =
            spinTypesVectorToJson(electronsVectorCollection.spinTypesVector())["SpinTypesVector"];
    j["ElectronsVectorCollection"]["PositionsVectorCollection"] =
            positionsVectorCollectionToJson(electronsVectorCollection.positionsVectorCollection())["PositionsVectorCollection"];
    return j;
}

nlohmann::json CollectionParser::atomsAndElectronsVectorToJson(const AtomsVector &atomsVector,
                                                               const ElectronsVector &electronsVector) {
    nlohmann::json j;
    j["AtomsVector"] = atomsVectorToJson(atomsVector)["AtomsVector"];
    j["ElectronsVector"] = electronsVectorToJson(electronsVector)["ElectronsVector"];
    return j;
}

PositionsVector CollectionParser::positionsVectorFromJson(const nlohmann::json &json) const {

    assert(json["PositionsVector"].is_array() && "File must have a PositionsVector array.");

    return arrayToPositionsVector(json["PositionsVector"]);
}

PositionsVector CollectionParser::arrayToPositionsVector(const nlohmann::json &json) const {
    std::vector<double> v = json.get<std::vector<double>>();
    assert(v.size()%3 == 0);
    return PositionsVector{Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v.data(), v.size())};
}

ElementTypesVector CollectionParser::elementTypesVectorFromJson(const nlohmann::json &json) const {
    assert(json["ElementTypesVector"].is_array());

    auto elementSymbols = json["AtomsVector"]["ElementTypesVector"].get<std::vector<std::string>>();
    ElementTypesVector elementTypesVector;
    for (auto &elementSymbol : elementSymbols) {
        elementTypesVector.append(Elements::ElementInfo::elementTypeForSymbol(elementSymbol));
    }
    return elementTypesVector;
}

SpinTypesVector CollectionParser::spinTypesVectorFromJson(const nlohmann::json &json) const {
    assert(json["SpinTypesVector"].is_array());

    auto spins = json["SpinTypesVector"].get<std::vector<int>>();
    SpinTypesVector spinTypesVector;
    for (int &spin : spins) {
        spinTypesVector.append( (Spin::SpinType) spin);
    }
    return spinTypesVector;
}

AtomsVector CollectionParser::atomsVectorFromJson(const nlohmann::json &json) const {
    assert(json["AtomsVector"].is_object() && "File must have a AtomsVector object.");

    return AtomsVector(positionsVectorFromJson(json["AtomsVector"]),
                       elementTypesVectorFromJson(json["AtomsVector"]));
}

ElectronsVector CollectionParser::electronsVectorFromJson(const nlohmann::json &json) const {
    assert(json["ElectronsVector"].is_object() && "File must have a ElectronsVector object.");

    return ElectronsVector(positionsVectorFromJson(json["ElectronsVector"]),
                           spinTypesVectorFromJson(json["ElectronsVector"]));
}

PositionsVectorCollection CollectionParser::positionsVectorCollectionFromJson(const nlohmann::json &json) const {
    assert(json["PositionsVectorCollection"].is_array());

    auto arrayOfArrays = json["PositionsVectorCollection"].get<std::vector<std::vector<double>>>();
    SpinTypesVector spinTypesVector;

    PositionsVectorCollection positionsVectorCollection;
    for (auto &array : arrayOfArrays) {
        positionsVectorCollection.append(arrayToPositionsVector(array));
    }

    return positionsVectorCollection;
}

ElectronsVectorCollection CollectionParser::electronsVectorCollectionFromJson(const nlohmann::json &json) const {
    assert(json["ElectronsVectorCollection"].is_array() && "File must have a ElectronsVectorCollection array.");

    return ElectronsVectorCollection(
            positionsVectorCollectionFromJson(json["ElectronsVectorCollection"]),
            spinTypesVectorFromJson(json["ElectronsVectorCollection"]));
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
