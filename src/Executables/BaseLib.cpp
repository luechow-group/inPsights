//
// Created by Michael Heuer on 21.06.18.
//

#include <omp.h>
#include <ParticlesVectorCollection.h>
#include <Particle.h>
#include <ElementType.h>
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <TestMolecules.h>
#include <Serialization.h>

int main(int argc, char *argv[]) {

    auto mol = TestMolecules::HeH::ElectronsInCores::normal;
    AtomsVectorCollection avc;
    avc.append(mol.atoms());
    avc.append(mol.atoms());

    //TEST EXPORT/IMPORT
    auto out = Serialization::yamlStringFrom<MolecularGeometry>(mol);
    std::cout << out << std::endl;
    std::string filename = "res.yaml";
    std::ofstream file(filename);
    file << out;
    file.close();
    auto y = YAML::LoadFile(filename);
    std::cout << y["Atoms"].as<AtomsVector>() << std::endl;
    std::cout << y["Electrons"].as<ElectronsVector>() << std::endl;


    auto out1 = Serialization::yamlStringFrom<MolecularGeometry>("Molecule",mol);
    std::cout << out1 << std::endl;

    auto out2 = Serialization::jsonStringFrom<MolecularGeometry>("mol",mol);
    std::cout << out2 << std::endl;

    auto out3 =  Serialization::yamlStringFrom<AtomsVectorCollection>("avc",avc);
    std::cout << out3 << std::endl;

    auto out4 = Serialization::jsonStringFrom<AtomsVectorCollection>("avc",avc);
    std::cout << out4 << std::endl;

    YAML::Node node;
    node["scalar"] = 2;
    node["Atoms"] = avc[0];

    Eigen::Matrix2d mat;
    mat << 1,2,3,4;

    YAML::Emitter emitter;
    emitter << YAML::Flow
            << YAML::BeginMap
            << YAML::Key << "Data" << YAML::Value << node
            << YAML::Key << "Matrix" << YAML::Value << 1234//mat
            << YAML::EndMap;
    std::cout << emitter.c_str() << std::endl;
}