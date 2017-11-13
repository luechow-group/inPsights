//
// Created by Morian Sonnet on 16.05.2017.
//

#include "Particle.h"
#include <cmath>

Particle::Particle(double &x, double &y, double &z): position{x,y,z} {
}

Particle::~Particle() {
}

double Particle::distance(const Particle &p1, const Particle &p2) {
    return (p1.position-p2.position).norm();
}

const Eigen::Vector3d& Particle::getPosition() const {
    return this->position;
}

void Particle::setPosition(double &x, double &y, double &z){
    this->position << x,y,z;
}
