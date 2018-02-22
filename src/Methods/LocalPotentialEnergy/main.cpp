//
// Created by leonard on 21.02.18.
//

#include <iostream>
#include "Particle.h"

int main(){
    std::cout << "hello world" << std::endl;

    auto position = Eigen::Vector3d(0.1,0.2,0.35);
    auto particle = Particle(position);

    std::cout << particle << std::endl;

    return 0;
}
