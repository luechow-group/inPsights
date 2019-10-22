/* Copyright (C) 2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef INPSIGHTS_MOTIF_H
#define INPSIGHTS_MOTIF_H

#include <list>
#include <string>
#include <Eigen/Core>
#include <yaml-cpp/yaml.h>

enum MotifType {unassigned=0, Core, CoreValence, Valence};

std::string toString(MotifType type);
MotifType fromString(const std::string& string);

class Motif{
public:
    Motif() = default;
    Motif(const std::list<Eigen::Index>& electronIndices, MotifType type = MotifType::unassigned);
    Motif(const std::list<Eigen::Index>& electronIndices,
          const std::list<Eigen::Index>& atomIndices,
          MotifType type = MotifType::unassigned);

    bool containsQ(Eigen::Index i) const;

    // needed for maps
    bool operator<(const Motif &rhs) const;
    bool operator>(const Motif &rhs) const;
    bool operator<=(const Motif &rhs) const;
    bool operator>=(const Motif &rhs) const;

    MotifType type() const;
    void setType(MotifType type_);

    const std::list<Eigen::Index> &electronIndices() const;
    void setElectronIndices(const std::list<Eigen::Index> &electronIndices);

    const std::list<Eigen::Index> &atomIndices() const;
    void setAtomIndices(const std::list<Eigen::Index> &atomIndices);

private:
    MotifType type_;
    std::list<Eigen::Index> electronIndices_;
    std::list<Eigen::Index> atomIndices_;
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
