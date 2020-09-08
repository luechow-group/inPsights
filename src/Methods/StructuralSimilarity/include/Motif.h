// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_MOTIF_H
#define INPSIGHTS_MOTIF_H

#include <set>
#include <string>
#include <Eigen/Core>
#include <yaml-cpp/yaml.h>
#include <MolecularSelection.h>


enum MotifType {unassigned=0, Core, Valence, CoreValence};

std::string toString(MotifType type);
MotifType fromString(const std::string& string);

class Motif : public MolecularSelection{
public:
    Motif() = default;
    Motif(ParticleIndices  electronIndices, MotifType type = MotifType::unassigned);
    Motif(ParticleIndices  electronIndices,
          ParticleIndices  nucleiIndices,
          MotifType type = MotifType::unassigned);

    // needed for maps
    bool operator<(const Motif &rhs) const;
    bool operator>(const Motif &rhs) const;
    bool operator<=(const Motif &rhs) const;
    bool operator>=(const Motif &rhs) const;

    MotifType type() const;
    void setType(MotifType type_);

private:
    MotifType type_;
};

namespace YAML {
    class Node; class Emitter;
    template<> struct convert<Motif> {
        static Node encode(const Motif &rhs);
        static bool decode(const Node &node, Motif &rhs);
    };
    Emitter &operator<<(Emitter &out, const Motif &rhs) ;
}

#endif //INPSIGHTS_MOTIF_H
