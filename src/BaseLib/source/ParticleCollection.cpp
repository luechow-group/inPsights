//
// Created by Michael Heuer on 29.10.17.
//

#include "ParticleCollection.h"

using namespace Eigen;

ParticleCollection::ParticleCollection()
        : AbstractCollection(),
          positions_(0)
{}

ParticleCollection::ParticleCollection(const VectorXd &positions) {
    auto size = (unsigned long) positions.size();
    assert(size >= 0 && "Vector is empty");
    assert(size%3 == 0 && "Vector is not 3N-dimensional");

    setNumberOfEntietes(size/3);
    positions_ = positions;
}

unsigned long ParticleCollection::numberOfParticles() const {
    return numberOfEntities_;
}

long ParticleCollection::calculateStartIndex(long i) const {
    assert(i <= long(numberOfEntities_) && "index is out of bounds");
    assert(i >= -long(numberOfEntities_) && "reverse index is out of bounds");
    if (i >= 0) return i*3;
    return (long(numberOfEntities_) + i)*3;
}

Particle ParticleCollection::operator[](long i) const {
    assert(i < long(numberOfEntities_) && "index is out of bounds");
    auto start = calculateStartIndex(i);
    auto position = positions_.segment(start,3);
    return Particle(position);
}

void ParticleCollection::insert(const Particle &particle, long i) {
    long start = calculateStartIndex(i);

    VectorXd before = positions_.head(start);
    VectorXd after = positions_.tail(numberOfEntities_*3-start);

    positions_.resize(numberOfEntities_*3+3);

    positions_.head(start) = before;
    positions_.segment(start,3) = particle.position();
    positions_.tail(numberOfEntities_*3-start) = after;

    ++numberOfEntities_;
}

std::ostream& operator<<(std::ostream& os, const ParticleCollection& pc){
    auto vec = pc.positionsAsEigenVector();
    Eigen::Map<Eigen::Matrix3Xd> mat(vec.data(),3,pc.numberOfParticles());
    os << static_cast<AbstractCollection>(pc)
       << mat.format(ParticleFormat::particleFormat) << std::endl;
    return os;
}

void ParticleCollection::prepend(const Particle &particle) {
    this->insert(particle,0);
}

void ParticleCollection::append(const Particle &particle) {
    this->insert(particle,numberOfEntities_);
}

Eigen::VectorXd ParticleCollection::positionsAsEigenVector() const {
    return positions_;
}

void ParticleCollection::permute(long i, long j) {
    assert( i >= 0 && i < numberOfEntities_
            && "Index i must be greater than zero and smaller than the number of particles." );
    assert( j >= 0 && j < numberOfEntities_
            && "Index j must be greater than zero and smaller than the number of particles." );
    if(i != j) {
        Eigen::Vector3d temp = positions_.segment(i*3,3);
        positions_.segment(i*3,3) = positions_.segment(j*3,3);
        positions_.segment(j*3,3) = temp;
    }
}
