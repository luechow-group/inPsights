//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_PARTICLESVECTOR_H
#define AMOLQCPP_PARTICLESVECTOR_H

#include "AbstractVector.h"
#include "Particle.h"
#include "PositionsVector.h"
#include "TypesVector.h"
#include <vector>

template<typename Type>
class ParticlesVector : public AbstractVector{
public:

    ParticlesVector()
            : AbstractVector(0),
              positionsVector_(),
              typesVector_(0)
    {}

    ParticlesVector(const PositionsVector &positionsVector)
            : AbstractVector(positionsVector.numberOfEntities()),
              positionsVector_(positionsVector),
              typesVector_(numberOfEntities())
    {}

    ParticlesVector(const PositionsVector &positionsVector,
                    const TypesVector<Type> &typesVector)
            : AbstractVector(positionsVector.numberOfEntities()),
              positionsVector_(positionsVector),
              typesVector_(typesVector) {
        assert(numberOfEntities() == positionsVector_.numberOfEntities()
               && numberOfEntities() == typesVector_.numberOfEntities()
               && "The number of entities in ParticlesVector, PositionsVector, and TypesVector must match.");
    }

    ParticlesVector(std::vector<Particle<Type>> particles)
            : ParticlesVector() {
        for (const auto& particle : particles){
            append(particle);
        }
    }
    
    Particle<Type> operator[](long i) const {
        return {typesVector_[i],positionsVector_[i]};
    }

    const PositionsVector & positionsVector() const {
        return positionsVector_;
    }

    PositionsVector & positionsVector(){
        return positionsVector_;
    }

    const TypesVector<Type> & typesVector() const {
        return typesVector_;
    }

    TypesVector<Type> & typesVector(){
        return typesVector_;
    }

    void insert(const Particle<Type>& particle, long i) {
        assert(i >= 0 && "The index must be positive.");
        assert(i <= numberOfEntities() && "The index must be smaller than the number of entities.");

        positionsVector_.insert(particle.position(),i);
        typesVector_.insert(particle.type(),i);
        incrementNumberOfEntities();

        assert(positionsVector_.numberOfEntities() == numberOfEntities());
        assert(typesVector_.numberOfEntities() == numberOfEntities());
    }

    void prepend(const Particle<Type> & particle) {
        this->insert(particle,0);
    }

    void append(const Particle<Type> & particle) {
        this->insert(particle,numberOfEntities());
    }

    void permute(long i, long j) override {
        positionsVector_.permute(i,j);
        typesVector_.permute(i,j);
    }

    friend std::ostream& operator<<(std::ostream& os, const ParticlesVector<Type> & pv){
        for (unsigned long i = 0; i < pv.numberOfEntities(); i++) {
            os << ToString::unsignedLongToString(i + 1) << " " << pv[i] << std::endl;
        }
        return os;
    }

protected:
    PositionsVector positionsVector_;
    TypesVector<Type> typesVector_;
};

using ElectronsVector = ParticlesVector<Spin>;
using AtomsVector = ParticlesVector<Element>;


#endif //AMOLQCPP_PARTICLESVECTOR_H
