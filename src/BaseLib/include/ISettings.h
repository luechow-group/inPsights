//
// Created by Michael Heuer on 2018-12-24.
//

#ifndef INPSIGHTS_ISETTINGS_H
#define INPSIGHTS_ISETTINGS_H

#include <Property.h>
#include <Varname.h>

#define YAML_SETTINGS_DECLARATION(classename)                   \
namespace YAML {                                                \
    class Node; class Emitter;                                  \
    template <typename Type> struct convert;                    \
    template<> struct convert<classename> {                     \
        static Node encode(const classename &rhs);              \
        static bool decode(const Node &node, classename &rhs);  \
    };                                                          \
    Emitter &operator<<(Emitter &out, const classename &rhs) ;  \
}

#define YAML_SETTINGS_DEFINITION(classname)                                 \
namespace YAML {                                                            \
    Node convert<classname>::encode(const classname &rhs) {                 \
        YAML::Node node;                                                    \
        rhs.appendToNode(node);                                             \
        return node;                                                        \
    }                                                                       \
    bool convert<classname>::decode(const Node &node, classname &rhs) {     \
        rhs = classname(node);                                              \
        return true;                                                        \
    }                                                                       \
    Emitter &operator<<(Emitter &out, const classname &rhs) {               \
        out << convert<classname>::encode(rhs);                             \
        return out;                                                         \
    }                                                                       \
}

namespace YAML {
    class Node;
}

namespace Settings{
    class ISettings {
    public:
        ISettings(const std::string& className);

        std::string name() const;

    protected:
        std::string className;
        virtual void appendToNode(YAML::Node &node) const = 0;
    };
}

#endif //INPSIGHTS_ISETTINGS_H
