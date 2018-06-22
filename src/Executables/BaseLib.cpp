//
// Created by Michael Heuer on 21.06.18.
//

#include <omp.h>
#include <ParticlesVectorCollection.h>
#include <ElementType.h>
#include <yaml-cpp/yaml.h>
#include <fstream>

int main(int argc, char *argv[]) {
    Eigen::Vector3d vec(1,2,3);
    Element et = Element::Ba;

    TypedParticle p{int(et),vec};
    Atom  a{et,vec};
    Electron e{Spin::alpha,vec};

    AtomsVector av{{a,a}};
    AtomsVectorCollection avc;
    avc.append(av);
    avc.append(av);

    YAML::Node n;
    n["test"] = av;
    n["position"] = av.positionsVector();
    n["pvc"] = avc.positionsVectorCollection();
    n["avc"] = avc;

    std::cout << n << std::endl;

    YAML::Emitter out;
    out << YAML::DoubleQuoted;
    //out << YAML::Flow; //JSON
    out << YAML::BeginMap;
    out << YAML::Key << "Atoms" << YAML::Value << av;
    out << YAML::Key << "Positions" << YAML::Value << av.positionsVector();
    out << YAML::Key << "PositionsColl" << YAML::Value << avc.positionsVectorCollection();
    out << YAML::Key << "AtomsColl" << YAML::Value << avc;
    out << YAML::EndMap;


    std::cout << "out\n" << out.c_str() << std::endl;

    std::ofstream file("results.yaml");
    file << out.c_str() << std::endl;

    std::ifstream in("results.yaml");
    std::cout << "in\n"<< YAML::Load(in);

}