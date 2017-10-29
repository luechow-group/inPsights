//
// Created by Michael Heuer on 29.10.17.
//

#include "Particle.h"

Particle::Particle(const Vector3d & position)
        : position_(position)
{}

Vector3d Particle::position() {
    return position_;
}

void Particle::position(const Vector3d &position) {
    position_ = position;
}
