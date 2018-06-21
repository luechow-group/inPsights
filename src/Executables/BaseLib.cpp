//
// Created by Michael Heuer on 21.06.18.
//

#include <omp.h>
#include <ParticlesVector.h>
#include <yaml-cpp/yaml.h>

int main(int argc, char *argv[]) {
    Eigen::Vector3d vec(1,2,3);
    Element et = Element::Ba;

    TypedParticle p{int(et),vec};
    Atom  a{et,vec};
    Electron e{Spin::alpha,vec};

    AtomsVector av{{a,a}};


    YAML::Emitter out;
    //out << YAML::DoubleQuoted << YAML::Flow;
    //out << YAML::BeginMap;
    //out << YAML::Key << "clusters";
    out << av;
    //out << YAML::EndMap;
    std::cout << out.c_str() << std::endl;
}