//
// Created by Michael Heuer on 2019-02-06.
//

#include "MotifEnergies.h"

void MotifEnergies::addPair(const std::pair<Motif, Motif>& pair, double energy){
    interactionEnergies_.emplace(std::make_pair(pair, energy));
}

const std::map<std::pair<Motif, Motif>, double>& MotifEnergies::interactionEnergies() const{
    return interactionEnergies_;
}

std::map<std::pair<Motif, Motif>, double>& MotifEnergies::interactionEnergies(){
    return interactionEnergies_;
}
