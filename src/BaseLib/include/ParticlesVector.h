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
#include <yaml-cpp/yaml.h>

template<typename Type>
class ParticlesVector : public ISliceable{
public:

    ParticlesVector()
            : ISliceable(0),
              positionsVector_(),
              typesVector_(0)
    {}

    ParticlesVector(const PositionsVector &positionsVector)
            : ISliceable(positionsVector.numberOfEntities()),
              positionsVector_(positionsVector),
              typesVector_(numberOfEntities())
    {}

    ParticlesVector(const PositionsVector &positionsVector,
                    const TypesVector<Type> &typesVector)
            : ISliceable(positionsVector.numberOfEntities()),
              positionsVector_(positionsVector),
              typesVector_(typesVector) {
        assert(numberOfEntities() == positionsVector_.numberOfEntities()
               && numberOfEntities() == typesVector_.numberOfEntities()
               && "The number of entities in ParticlesVector, PositionsVector, and TypesVector must match.");
    }

    ParticlesVector& slice(const Interval& interval, const Reset& resetType = Reset::Automatic /*TODO makes no sense here*/) {
        setSlice(interval,resetType);
        typesVector().slice(interval,resetType);
        positionsVector().slice(interval,resetType);

        //TODO CAREFUL WITH RESET
        return *this;
    }

    ParticlesVector& entity(long i, const Reset& resetType = Reset::Automatic) {
        return slice(Interval(i), resetType);
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

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation, const Usage &usage) {
        positionsVector_.permute(permutation,usage);
        typesVector_.permute(permutation,usage);
    }

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) override {
        positionsVector_.permute(permutation);
        typesVector_.permute(permutation);
    }


    friend std::ostream& operator<<(std::ostream& os, const ParticlesVector<Type> & pv){
        //TODO print only slice?
        for (long i = 0; i < pv.numberOfEntities(); i++) {
            os << ToString::longToString(i + 1) << " " << pv[i] << std::endl;
        }
        return os;
    }

    bool operator==(const ParticlesVector<Type> &other) const {
        return ISliceable::operator==(other)
                && (typesVector() == other.typesVector())
                && (positionsVector() == other.positionsVector());
    }

    bool operator!=(const ParticlesVector<Type> &other) const {
        return !(*this == other);
    }



protected:
    PositionsVector positionsVector_;
    TypesVector<Type> typesVector_;
};

using TypedParticlesVector = ParticlesVector<int>;
using ElectronsVector = ParticlesVector<Spin>;
using AtomsVector = ParticlesVector<Element>;

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

#endif //AMOLQCPP_PARTICLESVECTOR_H
