//
// Created by Michael Heuer on 29.10.17.
//

#include "ParticleCollection.h"

ParticleCollection::ParticleCollection()
        : numberOfParticles_(0),
          positions_(0)
{}

ParticleCollection::ParticleCollection(const VectorXd &positions) {
    long size = positions.size();
    assert(size >= 0 && "Vector is empty");
    assert(size%3 == 0 && "Vector is not 3N-dimensional");

    numberOfParticles_ = size/3;
    positions_ = positions;
}

long ParticleCollection::numberOfParticles() {
    return numberOfParticles_;
}

long ParticleCollection::calculateStartIndex(long i) {
    assert(i <= numberOfParticles_ && "index is out of bounds");
    assert(i >= -numberOfParticles_ && "reverse index is out of bounds");
    if (i >= 0) return i*3;
    return (numberOfParticles_ + i)*3;
}

Particle ParticleCollection::operator[](long i) {
    assert(i < numberOfParticles_ && "index is out of bounds");
    long start = calculateStartIndex(i);
    Vector3d position = positions_.segment(start,3);
    return Particle(position);
}

void ParticleCollection::insert(const Particle &particle, long i) {
    long start = calculateStartIndex(i);

    VectorXd before = positions_.head(start);
    VectorXd after = positions_.tail(numberOfParticles_*3-start);

    positions_.resize(numberOfParticles_*3+3);
    positions_ << before, particle.position(), after;
    ++numberOfParticles_;
}

void ParticleCollection::prepend(const Particle &particle) {
    this->insert(particle,0);
}

void ParticleCollection::append(const Particle &particle) {
    this->insert(particle,numberOfParticles_);
}

Eigen::VectorXd ParticleCollection::positionsAsEigenVector() {
    return positions_;
}
