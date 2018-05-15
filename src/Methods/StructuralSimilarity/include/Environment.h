//
// Created by Michael Heuer on 03.05.18.
//

#ifndef AMOLQCPP_ENVIRONMENT_H
#define AMOLQCPP_ENVIRONMENT_H

#include <Eigen/Core>
#include <ParticlesVector.h>

class Environment{
public:
    Environment(AtomsVector atoms, unsigned centerId)
            : atoms_(atoms),
              centerId_(centerId){}

    AtomsVector atoms_;
    unsigned centerId_;
};

#endif //AMOLQCPP_ENVIRONMENT_H
