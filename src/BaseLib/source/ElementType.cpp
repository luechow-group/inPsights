//
// Created by Michael Heuer on 05.02.18.
//

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
