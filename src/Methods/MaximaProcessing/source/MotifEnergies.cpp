//
// Created by Michael Heuer on 2019-02-06.
//

#include "MotifEnergies.h"
#include <iostream>
void MotifEnergies::addPair(const std::vector<long>& motifIds, double energy){
    interactionEnergies_.emplace(std::make_pair(motifIds, energy));
}

const std::map<std::vector<long>, double>& MotifEnergies::interactionEnergies() const{
    return interactionEnergies_;
}

std::map<std::vector<long>, double>& MotifEnergies::interactionEnergies(){
    return interactionEnergies_;
}

namespace YAML {
    Node convert<MotifEnergies>::encode(const MotifEnergies &rhs) {
        Node mainNode;

        for(const auto& motifEnergyPair : rhs.interactionEnergies()) {
            Node node;
            node["MotifIds"] = motifEnergyPair.first;
            node["Energy"] = motifEnergyPair.second;

            mainNode.push_back(node);
        }

        return mainNode;
    }

    bool convert<MotifEnergies>::decode(const Node &node, MotifEnergies &rhs) {
        if (!node.IsSequence())
            return false;

        MotifEnergies motifEnergies;

        for (const auto &n : node)
            motifEnergies.addPair(
                    n["MotifIds"].as<std::vector<long>>(),
                    n["Energy"].as<double>());

        rhs = motifEnergies;

        return true;
    }

    Emitter &operator<<(Emitter &out, const MotifEnergies &rhs) {
        out << BeginSeq;

        for(const auto& i : rhs.interactionEnergies()){
            out << BeginMap
            << Key << "MotifIds" << Flow << Value << i.first
            << Key << "Energy" <<  Value << i.second
            << EndMap;
        }

        out << EndSeq;
        return out;
    };
}
