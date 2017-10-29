//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_PARTICLE_H
#define AMOLQCGUI_PARTICLE_H

#include <Eigen/Core>

using namespace Eigen;

class Particle {
public:
    explicit Particle(const Vector3d& position);

    Vector3d position();
    void position(const Vector3d& position);

protected:
    Vector3d position_;
};
#endif //AMOLQCGUI_PARTICLE_H
