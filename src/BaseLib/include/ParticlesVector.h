//
// Created by Michael Heuer on 29.10.17.
//

#ifndef INPSIGHTS_PARTICLESVECTOR_H
#define INPSIGHTS_PARTICLESVECTOR_H

#include "AbstractVector.h"
#include "Particle.h"
#include "PositionsVector.h"
#include "TypesVector.h"
#include <yaml-cpp/yaml.h>
#include <vector>
#include <memory>

template<typename Type>
class ParticlesVector : public AbstractVector{
public:

    ParticlesVector()
            : AbstractVector(0),
              positionsVector_(),
              typesVector_(0),
              linkedParticles_(0) {
        initializeLinkedParticles();
    }

    ParticlesVector(const PositionsVector &positionsVector)
            : AbstractVector(positionsVector.numberOfEntities()),
              positionsVector_(positionsVector),
              typesVector_(numberOfEntities()),
              linkedParticles_(0) {
        initializeLinkedParticles();
    }

    ParticlesVector(const PositionsVector &positionsVector,
                    const TypesVector<Type> &typesVector)
            : AbstractVector(positionsVector.numberOfEntities()),
              positionsVector_(positionsVector),
              typesVector_(typesVector),
              linkedParticles_(0) {
        assert(numberOfEntities() == positionsVector_.numberOfEntities()
               && numberOfEntities() == typesVector_.numberOfEntities()
               && "The number of entities in ParticlesVector, PositionsVector, and TypesVector must match.");
        initializeLinkedParticles();
    }

    ParticlesVector(const ParticlesVector& pv)
            : ParticlesVector(pv.positionsVector(),pv.typesVector()) {}

    void initializeLinkedParticles(){
        for (long i = 0; i < numberOfEntities(); ++i) {
            linkedParticles_.emplace_back(std::make_shared<LinkedParticle<Type>>(
                    positionsVector_.dataRef(i), &typesVector_.dataRef(i)[0]));//TODO or give LinkedParticle a VectorXi(1) ref?
        }
    }

    ParticlesVector(std::vector<Particle<Type>> particles)
            : ParticlesVector() {
        for (const auto& particle : particles){
            append(particle);
        }
        initializeLinkedParticles();
    }

    std::shared_ptr<LinkedParticle<Type>> linkedParticle(long i) const {
        return linkedParticles_[i];
    }

    Particle<Type> operator[](long i) const {
        return {typesVector_[i],positionsVector_[i]};
    }

    ParticlesVector<Type> operator[](std::list<long> indices) const {
        ParticlesVector newVector;
        for (long index : indices) {
            newVector.append(Particle<Type>{typesVector_[index],positionsVector_[index]});
        }
        return newVector;
    }

    ParticlesVector<Type> getFirstElements(long i) {
        assert(i >= 0);

        ParticlesVector newVector(
                PositionsVector(positionsVector_.asEigenVector().head(3*i)),
                TypesVector<Type>(typesVector_.asEigenVector().head(i))
        );
        return newVector;
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

        linkedParticles_.resize(0); //TODO find alternative to recreate the entire vector (what happens to existing linked particles)
        initializeLinkedParticles();
        assert(long(linkedParticles_.size()) == numberOfEntities());
    }

    void prepend(const Particle<Type> & particle) {
        this->insert(particle,0);
    }

    void append(const Particle<Type> & particle) {
        this->insert(particle,numberOfEntities());
    }

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) override {
        positionsVector_.permute(permutation);
        typesVector_.permute(permutation);
    }


    friend std::ostream& operator<<(std::ostream& os, const ParticlesVector<Type> & pv){
        for (long i = 0; i < pv.numberOfEntities(); i++) {
            os << ToString::longToString(i + 1) << " " << pv[i] << std::endl;
        }
        return os;
    }

    bool operator==(const ParticlesVector<Type> &other) const {
        return (typesVector() == other.typesVector())
        && (positionsVector() == other.positionsVector());
    }

    bool operator!=(const ParticlesVector<Type> &other) const {
        return !(*this == other);
    }

protected:
    PositionsVector positionsVector_;
    TypesVector<Type> typesVector_;
public:
    std::vector<std::shared_ptr<LinkedParticle<Type>>> linkedParticles_;
};

using TypedParticlesVector = ParticlesVector<int>;
using ElectronsVector = ParticlesVector<Spin>;
using AtomsVector = ParticlesVector<Element>;

using LinkedTypedParticle = LinkedParticle<int>;
using LinkedElectron = LinkedParticle<Spin>;
using LinkedAtoms = LinkedParticle<Element>;

namespace YAML {
    template<typename Type> struct convert<ParticlesVector<Type>> {
        static Node encode(const ParticlesVector<Type> & pv){
            Node node;
            node["Types"] = pv.typesVector();
            node["Positions"] = pv.positionsVector();
            return node;

        }
        static bool decode(const Node& nodes, ParticlesVector<Type> & rhs){
            if(!nodes.IsMap())
                return false;
            ParticlesVector<Type> pv(
                    nodes["Positions"].as<PositionsVector>(),
                    nodes["Types"].as<TypesVector<Type>>());
            rhs = pv;
            return true;
        }
    };

    template<typename Type>
    Emitter& operator<< (Emitter& out, const ParticlesVector<Type>& pv){
        out << BeginMap// << Newline
            << Key << "Types" << Value << pv.typesVector() << Newline
            << Key << "Positions" << Value <<  pv.positionsVector()// << Newline
            << EndMap;
        return out;
    };
}

#endif //INPSIGHTS_PARTICLESVECTOR_H
