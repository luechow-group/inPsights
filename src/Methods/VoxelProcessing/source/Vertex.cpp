//
// Created by Michael Heuer on 2019-01-03.
//

#include <Vertex.h>
#include <yaml-cpp/yaml.h>
#include <EigenYamlConversion.h>

 Vertex::Vertex(const Eigen::Vector3f& position)
 : position(position),normal(Eigen::Vector3f::Zero()) {}

namespace YAML {
    Node convert<Vertex>::encode(const Vertex &rhs) {
        Node node;
        node.push_back(rhs.position);
        node.push_back(rhs.normal);
        return node;
    }
    bool convert<Vertex>::decode(const Node &node, Vertex &rhs) {
     if (!node.IsSequence())
         return false;

     rhs.position = node[0].as<Eigen::Vector3f>();
     rhs.normal = node[1].as<Eigen::Vector3f>();
     return true;
    }

    Emitter &operator<<(Emitter &out, const Vertex &rhs) {
        out << Flow << BeginSeq << rhs.position << rhs.normal << EndSeq;
        return out;
    }
}
