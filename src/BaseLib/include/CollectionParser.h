//
// Created by Michael Heuer on 08.12.17.
//

#ifndef AMOLQCPP_COLLECTIONPARSER_H
#define AMOLQCPP_COLLECTIONPARSER_H

#include "ElectronsVectorCollection.h"
#include "AtomsVector.h"

#include <nlohmann/json.hpp>

class CollectionParser{
public:

    nlohmann::json positionsVectorToJson(const PositionsVector &positionsVector);
    nlohmann::json elementTypesVectorToJson(const ElementTypesVector& elementTypesVector);
    nlohmann::json spinTypesVectorToJson(const SpinTypesVector& spinTypesVector);

    nlohmann::json atomsVectorToJson(const AtomsVector& atomsVector);
    nlohmann::json electronsVectorToJson(const ElectronsVector& electronsVector);

    nlohmann::json positionsVectorCollectionToJson(const PositionsVectorCollection& positionsVectorCollection);
    nlohmann::json electronsVectorCollectionToJson(const ElectronsVectorCollection& electronsVectorCollection);

    nlohmann::json atomsAndElectronsVectorToJson(const AtomsVector& atomsVector,const ElectronsVector& electronsVector);


    PositionsVector positionsVectorFromJson(const nlohmann::json& json) const;
    ElementTypesVector elementTypesVectorFromJson(const nlohmann::json& json) const;
    SpinTypesVector spinTypesVectorFromJson(const nlohmann::json& json) const;

    AtomsVector atomsVectorFromJson(const nlohmann::json& json) const;
    ElectronsVector electronsVectorFromJson(const nlohmann::json& json) const;

    PositionsVectorCollection positionsVectorCollectionFromJson(const nlohmann::json& json) const;
    ElectronsVectorCollection electronsVectorCollectionFromJson(const nlohmann::json& json) const;

    nlohmann::json json, readJSON(const std::string& filename);
    void writeJSON(const nlohmann::json& json, const std::string& filename);

private:
    PositionsVector arrayToPositionsVector(const nlohmann::json &json) const;
};

#endif //AMOLQCPP_COLLECTIONPARSER_H
