//
// Created by Michael Heuer on 29.10.17.
//

#include "Particle.h"

using namespace Eigen;

Particle::Particle(const Vector3d & position)
        : position_(position) {}

Particle::Particle(double x, double y, double z)
        : position_(x,y,z) {}

Vector3d Particle::position() const {
    return position_;
}

void Particle::position(const Vector3d &position) {
    position_ = position;
}

double Particle::distance(const Particle &p1, const Particle &p2) {
    return (p1.position()-p2.position()).norm();
}

std::ostream& operator<< (std::ostream& os, const Particle& p) {
    os << p.position().format(ParticleFormat::particleFormat) << std::endl;
    return os;
}
