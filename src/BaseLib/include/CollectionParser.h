//
// Created by Michael Heuer on 08.12.17.
//

#ifndef AMOLQCPP_COLLECTIONPARSER_H
#define AMOLQCPP_COLLECTIONPARSER_H

#include "ElectronCollection.h"
#include "AtomCollection.h"

#include <json.hpp>

//namespace CollectionParser{
class CollectionParser{
public:

    nlohmann::json atomCollectionToJson(const AtomCollection& atomCollection);
    nlohmann::json electronCollectionToJson(const ElectronCollection& electronCollection);

    nlohmann::json particleCollectionToJson(const ParticleCollection& particleCollection);
    void writeJSON(const nlohmann::json& json, const std::string& filename);
    //static ParticleCollection readParticleCollection(std::string & filename);
};

#endif //AMOLQCPP_COLLECTIONPARSER_H
