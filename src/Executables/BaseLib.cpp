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
    PositionsVectorCollection pvc;
    pvc.append(mol.atoms().positionsVector());
    pvc.append(mol.atoms().positionsVector());

    auto out = Serialization::yamlStringFrom<PositionsVectorCollection>("Pvc",pvc);
    auto out2 = Serialization::yamlStringFrom<MolecularGeometry>("Pvc",mol);


    std::cout << out2 << std::endl;

    std::ofstream filejson("results.json");
    filejson << out.c_str();
    filejson.close();
}