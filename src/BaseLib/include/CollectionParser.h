//
// Created by Michael Heuer on 08.12.17.
//

#ifndef AMOLQCPP_COLLECTIONPARSER_H
#define AMOLQCPP_COLLECTIONPARSER_H

#include "ElectronsVector.h"
#include "AtomsVector.h"

#include <nlohmann/json.hpp>

class CollectionParser{
public:

    nlohmann::json atomsVectorToJson(const AtomsVector& atomsVector);
    nlohmann::json electronsVectorToJson(const ElectronsVector& electronsVector);
    nlohmann::json positionsVectorToJson(const PositionsVector &positionsVector);

    AtomsVector atomsVectorFromJson (const std::string& filename);
    ElectronsVector electronsVectorFromJson(const std::string& filename);
    PositionsVector positionsVectorFromJson(const std::string &filename);

    nlohmann::json json, readJSON(const std::string& filename);
    void writeJSON(const nlohmann::json& json, const std::string& filename);

private:
    PositionsVector array2DToPositionsVector(const nlohmann::json &coordinates);
};

#endif //AMOLQCPP_COLLECTIONPARSER_H
