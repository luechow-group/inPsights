//
// Created by Michael Heuer on 21.06.18.
//

#include <omp.h>
#include <Particle.h>
#include <yaml-cpp/yaml.h>

int main(int argc, char *argv[]) {
    Eigen::Vector3d vec(1,2,3);
    Element et = Element::Ba;

    TypedParticle p{int(et),vec};
    Atom  a{et,vec};
    Electron e{Spin::alpha,vec};



    YAML::Emitter out;
    //out << YAML::DoubleQuoted << YAML::Flow;
    //out << YAML::BeginMap;
    //out << YAML::Key << "clusters";
    out << p;
    out << et;
    out << vec;
    out << a;
    out << e;
    //out << YAML::EndMap;
    std::cout << out.c_str() << std::endl;
}