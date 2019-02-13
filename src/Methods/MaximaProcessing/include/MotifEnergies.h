//
// Created by Michael Heuer on 2019-02-06.
//

#ifndef INPSIGHTS_MOTIFENERGIES_H
#define INPSIGHTS_MOTIFENERGIES_H

#include <Motif.h>
#include <map>

class MotifEnergies{
public:
    void addPair(const std::vector<long>& motifIds, double energy);

    const std::map<std::vector<long>, double>& interactionEnergies() const;

    std::map<std::vector<long>, double>& interactionEnergies();

private:
    std::map<std::vector<long>, double> interactionEnergies_;
};

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<>
    struct convert<MotifEnergies> {
        static Node encode(const MotifEnergies& rhs);
        static bool decode(const Node& node, MotifEnergies& rhs);
    };
    Emitter& operator<< (Emitter& out, const MotifEnergies& p);
}

#endif //INPSIGHTS_MOTIFENERGIES_H
