//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_PARTICLECOLLECTION_H
#define AMOLQCGUI_PARTICLECOLLECTION_H

#include "Particle.h"

using namespace Eigen;

class ParticleCollection{
public:

    explicit ParticleCollection(const VectorXd& positions);

    Particle operator[](long i);
    long size();

protected:
    long size_;
    VectorXd positions_;

private:
    long calculateStartIndex(long i);

};

#endif //AMOLQCGUI_PARTICLECOLLECTION_H
