//
// Created by Michael Heuer on 22.06.18.
//

#include <gtest/gtest.h>
#include <sstream>
#include <TestMolecules.h>
#include <nlohmann/json.hpp>
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <sstream>

using namespace testing;

class AYamlTest : public Test {
public:
};

TEST_F(AYamlTest, JsonCompatibility){

    MolecularGeometry mol = TestMolecules::HeH::ElectronsInCores::normal;

    YAML::Emitter out;
    out << YAML::Flow;
    out << YAML::DoubleQuoted;
    out << YAML::BeginMap;
    out << YAML::Key << "Molecule" << YAML::Value << mol;
    out << YAML::EndMap;

    std::ofstream ofstream("AYamlTest.json");

    std::string outstring = out.c_str();;
    ofstream << outstring;
    ofstream.close();

    std::cout << outstring<< std::endl;
    
    
    std::string expected = "{\"Molecule\": {\"Atoms\": [[\"He\", [0, 0, 0.37]], [\"H\", [0, 0, -0.37]]], \"Electrons\": [[\"a\", [0, 0, 0.37]], [\"a\", [0, 0, 0.37]], [\"b\", [0, 0, -0.37]]]}}";

    ASSERT_EQ(outstring,expected);
}