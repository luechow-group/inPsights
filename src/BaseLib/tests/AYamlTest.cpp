// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <gmock/gmock.h>
#include <sstream>
#include <TestMolecules.h>
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <sstream>

using namespace testing;

TEST(AYamlTest, JsonCompatibility){

    MolecularGeometry mol = TestMolecules::HeH::ElectronsInCores::normal;

    YAML::Emitter out;
    out << YAML::Flow;
    out << YAML::DoubleQuoted;
    out << YAML::BeginMap;
    out << YAML::Key << "Molecule" << YAML::Value << mol;
    out << YAML::EndMap;

    std::string outstring = out.c_str();
    std::string expected ="{\"Molecule\": {\"Atoms\": {\"Types\": [\"He\", \"H\"],\n  \"Positions\": [\n      [0, 0, 0.37],\n      [0, 0, -0.37],\n      ]}\n\"Electrons\": {\"Types\": [\"a\", \"a\", \"b\"],\n  \"Positions\": [\n      [0, 0, 0.37],\n      [0, 0, 0.37],\n      [0, 0, -0.37],\n      ]}}}";
    ASSERT_EQ(outstring,expected);
}