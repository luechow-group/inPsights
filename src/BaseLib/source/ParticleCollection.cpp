//
// Created by Michael Heuer on 29.10.17.
//

#include "ParticleCollection.h"

ParticleCollection::ParticleCollection(const VectorXd &positions) {
    long size = positions.size();
    assert(size >= 0 && "Vector is empty");
    assert(size%3 == 0 && "Vector is not 3N-dimensional");

    size_ = size/3;
    positions_ = positions;
}

long ParticleCollection::size() {
    return size_;
}

long ParticleCollection::calculateStartIndex(long i) {
    assert(i < size_ && "index is out of bounds");
    assert(i >= -size_ && "reverse index is out of bounds");
    if (i >= 0) return i*3;
    return (size_ + i)*3;
}

Particle ParticleCollection::operator[](long i) {
    long start = calculateStartIndex(i);
    Vector3d position = positions_.segment(start,3);
    return Particle(position);
}
