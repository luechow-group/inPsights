/* Copyright (C) 2018-2019 Michael Heuer.
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

#include "ElementType.h"
#include "ElementInfo.h"
#include <yaml-cpp/yaml.h>

Element Elements::first() {
    return Element::H;
};

Element Elements::last(){
    return Element::Og;
};

Element Elements::elementFromInt(int type){
    return static_cast<Element>(type);
};

int Elements::elementToInt(Element element){
    return int(element);
};

std::ostream& operator<< (std::ostream& os, const Element & e){
    os <<  Elements::ElementInfo::symbol(e);
    return os;
};

namespace YAML {
    Node convert<Element>::encode(const Element &rhs) {
        Node node = YAML::convert<std::string>::encode(Elements::ElementInfo::symbol(rhs));
        return node;
    }

    bool convert<Element>::decode(const Node &node, Element &rhs) {
        if (!node.IsScalar()) {
            return false;
        }
        rhs = Elements::ElementInfo::elementTypeFromSymbol(node.as<std::string>());
        return true;
    }

    Emitter &operator<<(Emitter &out, const Element &e) {
        out << Elements::ElementInfo::symbol(e);
        return out;
    }
}
