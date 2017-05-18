//
// Created by Moria on 16.05.2017.
//

#include "Particle.h"
#include <cmath>

Particle::Particle(double x, double y, double z): position{x,y,z} {
}

Particle::~Particle() {
}

double Particle::distance(const Particle &p1, const Particle &p2) {
    return sqrt(pow(p1.get_position()[0]-p2.get_position()[0],2)+pow(p1.get_position()[1]-p2.get_position()[1],2)+pow(p1.get_position()[2]-p2.get_position()[2],2));
}

const std::vector<double>& Particle::get_position() const {
    return this->position;
}
