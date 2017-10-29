//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_PARTICLECOLLECTION_H
#define AMOLQCGUI_PARTICLECOLLECTION_H

#include "Particle.h"
#include <vector>

using namespace Eigen;

class ParticleCollection{
public:

    explicit ParticleCollection(const VectorXd& positions);

    Particle operator[](long i);
    long size();


    void insert(const Particle& particle, long i);
    void append(const Particle& particle);
    void prepend(const Particle& particle);
    /* TODO
    void replace(long i);
    void remove(long i);
    ParticleCollection part(std::vector<long> indices);
    */

protected:
    long size_;
    VectorXd positions_;

private:
    long calculateStartIndex(long i);

};

#endif //AMOLQCGUI_PARTICLECOLLECTION_H
