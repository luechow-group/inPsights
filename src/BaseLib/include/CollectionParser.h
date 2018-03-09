//
// Created by Michael Heuer on 08.12.17.
//

#ifndef AMOLQCPP_COLLECTIONPARSER_H
#define AMOLQCPP_COLLECTIONPARSER_H

#include "ElectronCollection.h"
#include "AtomCollection.h"

#include <nlohmann/json.hpp>

class CollectionParser{
public:

    nlohmann::json atomCollectionToJson(const AtomCollection& atomCollection);
    nlohmann::json electronCollectionToJson(const ElectronCollection& electronCollection);
    nlohmann::json positionsVectorToJson(const PositionsVector &positionsVector);

    AtomCollection atomCollectionFromJson (const std::string& filename);
    ElectronCollection electronCollectionFromJson(const std::string& filename);
    PositionsVector positionsVectorFromJson(const std::string &filename);

    nlohmann::json json, readJSON(const std::string& filename);
    void writeJSON(const nlohmann::json& json, const std::string& filename);

private:
    PositionsVector array2DToPositionsVector(const nlohmann::json &coordinates);
};

#endif //AMOLQCPP_COLLECTIONPARSER_H
