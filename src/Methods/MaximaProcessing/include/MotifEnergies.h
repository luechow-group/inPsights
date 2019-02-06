//
// Created by Michael Heuer on 2019-02-06.
//

#ifndef INPSIGHTS_MOTIFENERGIES_H
#define INPSIGHTS_MOTIFENERGIES_H

#include <Motif.h>
#include <map>

class MotifEnergies{
public:
    void addPair(const std::pair<Motif, Motif>& pair, double energy);

    const std::map<std::pair<Motif, Motif>, double>& interactionEnergies() const;

    std::map<std::pair<Motif, Motif>, double>& interactionEnergies();

private:
    std::map<std::pair<Motif, Motif>, double> interactionEnergies_;
};


#endif //INPSIGHTS_MOTIFENERGIES_H
