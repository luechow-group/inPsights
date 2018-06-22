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
#include <nlohmann/json.hpp>
#include <CollectionParser.h>

int main(int argc, char *argv[]) {

    auto mol = TestMolecules::HeH::ElectronsInCores::normal;
    AtomsVectorCollection avc;
    avc.append(mol.atoms());
    avc.append(mol.atoms());

    auto out1 = Serialization::yamlStringFrom<MolecularGeometry>("mol",mol);
    std::cout << out1 << std::endl;

    auto out2 = Serialization::jsonStringFrom<MolecularGeometry>("mol",mol);
    std::cout << out2 << std::endl;

    auto out3 =  Serialization::yamlStringFrom<AtomsVectorCollection>("avc",avc);
    std::cout << out3 << std::endl;

    auto out4 = Serialization::jsonStringFrom<AtomsVectorCollection>("avc",avc);
    std::cout << out4 << std::endl;

}