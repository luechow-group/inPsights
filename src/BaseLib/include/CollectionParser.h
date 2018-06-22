//
// Created by Michael Heuer on 08.12.17.
//

#ifndef AMOLQCPP_COLLECTIONPARSER_H
#define AMOLQCPP_COLLECTIONPARSER_H

#include <Eigen/Eigenvalues>
#include <nlohmann/json.hpp>

#include "ParticlesVectorCollection.h"
#include "PositionsVectorCollection.h"
#include "fstream"
#include <yaml-cpp/yaml.h>

namespace Serialization{
    template<typename Type>
    std::string yamlStringFrom(const std::string &key, const Type &value){
        YAML::Emitter out;
        out << YAML::BeginMap;
        out << YAML::Key << key << YAML::Value << value;
        out << YAML::EndMap;
        return out.c_str();
    };

    template<typename Type>
    std::string yamlStringFrom(const Type &scalar){
        YAML::Emitter out;
        out << scalar;
        return out.c_str();
    };

    template<typename Type>
    std::string jsonStringFrom(const std::string &key, const Type &value){
        YAML::Emitter out;
        out << YAML::Flow;
        out << YAML::DoubleQuoted;
        out << YAML::BeginMap;
        out << YAML::Key << key << YAML::Value << value;
        out << YAML::EndMap;
        return out.c_str();
    };

    template<typename Type>
    std::string jsonStringFrom(const Type &scalar){
        YAML::Emitter out;
        out << YAML::Flow;
        out << YAML::DoubleQuoted;
        out << scalar;
        return out.c_str();
    };
};


namespace CollectionParser{
    nlohmann::json positionsVectorToJson(const PositionsVector &positionsVector);
    nlohmann::json elementTypesVectorToJson(const ElementTypesVector& elementTypesVector);
    nlohmann::json spinTypesVectorToJson(const SpinTypesVector& spinTypesVector);

    nlohmann::json atomsVectorToJson(const AtomsVector& atomsVector);
    nlohmann::json electronsVectorToJson(const ElectronsVector& electronsVector);

    nlohmann::json positionsVectorCollectionToJson(const PositionsVectorCollection& positionsVectorCollection);
    nlohmann::json electronsVectorCollectionToJson(const ElectronsVectorCollection& electronsVectorCollection);

    nlohmann::json atomsAndElectronsVectorToJson(const AtomsVector& atomsVector,const ElectronsVector& electronsVector);

    nlohmann::json selfAdjointEigenSolverResultsToJsonArray(const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> &eigenSolver);

    PositionsVector positionsVectorFromJson(const nlohmann::json& json);
    ElementTypesVector elementTypesVectorFromJson(const nlohmann::json& json);
    SpinTypesVector spinTypesVectorFromJson(const nlohmann::json& json) ;

    AtomsVector atomsVectorFromJson(const nlohmann::json& json) ;
    ElectronsVector electronsVectorFromJson(const nlohmann::json& json) ;

    PositionsVectorCollection positionsVectorCollectionFromJson(const nlohmann::json& json) ;
    ElectronsVectorCollection electronsVectorCollectionFromJson(const nlohmann::json& json) ;

    nlohmann::json readJSON(const std::string& filename);
    void writeJSON(const nlohmann::json& json, const std::string& filename);

    PositionsVector arrayToPositionsVector(const nlohmann::json &json) ;
};


#endif //AMOLQCPP_COLLECTIONPARSER_H
