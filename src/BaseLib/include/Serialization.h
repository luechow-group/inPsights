// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_COLLECTIONPARSER_H
#define INPSIGHTS_COLLECTIONPARSER_H

#include <yaml-cpp/yaml.h>
#include <fstream>

#include "ParticlesVectorCollection.h"
#include "PositionsVectorCollection.h"

namespace Serialization{

    template<typename Type>
    std::string yamlStringFrom(const std::string &key, const Type &value,
                               YAML::EMITTER_MANIP format = YAML::EMITTER_MANIP::Block){
        YAML::Emitter out;
        out << YAML::BeginMap << YAML::Newline;
        out << format;
        out << YAML::Key << key << YAML::Value << value << YAML::Newline;
        out << YAML::EndMap;
        return out.c_str();
    };

    template<typename Type>
    std::string yamlStringFrom(const Type &scalar,
                               YAML::EMITTER_MANIP format = YAML::EMITTER_MANIP::Block){
        YAML::Emitter out;
        out << format;
        out << scalar << YAML::Newline;
        return out.c_str();
    };

    template<typename Type>
    std::string jsonStringFrom(const std::string &key, const Type &value){
        YAML::Emitter out;
        out << YAML::Flow;
        out << YAML::DoubleQuoted;
        out << YAML::BeginMap;
        out << YAML::Key << key  << YAML::Newline << YAML::Value << value << YAML::Newline;
        out << YAML::EndMap;
        return out.c_str();
    };

    template<typename Type>
    std::string jsonStringFrom(const Type &scalar){
        YAML::Emitter out;
        out << YAML::Flow;
        out << YAML::DoubleQuoted;
        out << scalar << YAML::Newline;
        return out.c_str();
    };
};


#endif //INPSIGHTS_COLLECTIONPARSER_H
