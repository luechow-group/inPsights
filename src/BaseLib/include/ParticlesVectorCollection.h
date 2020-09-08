// Copyright (C) 2017-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_PARTICLESVECTORCOLLECTION_H
#define INPSIGHTS_PARTICLESVECTORCOLLECTION_H

#include <vector>
#include "ParticlesVector.h"
#include "PositionsVectorCollection.h"

template <typename Type>
class ParticlesVectorCollection : public AbstractVector{
public:
    ParticlesVectorCollection()
            : AbstractVector(0),
              positionsVectorCollection_(),
              typesVector_(0) {}

    ParticlesVectorCollection(const PositionsVectorCollection& positionsVectorCollection)
            : AbstractVector(positionsVectorCollection.numberOfEntities()),
              positionsVectorCollection_(positionsVectorCollection),
              typesVector_(positionsVectorCollection.numberOfPositionsEntities()) {}

    ParticlesVectorCollection(const TypesVector<Type> &typesVector)
            : AbstractVector(0),
              positionsVectorCollection_(),
              typesVector_(typesVector) {}

    ParticlesVectorCollection(const ParticlesVector<Type> &particlesVector)
            : ParticlesVectorCollection(std::vector<ParticlesVector<Type>>({particlesVector})){}

    ParticlesVectorCollection(const std::vector<ParticlesVector<Type>> &particlesVectorVector)
            : AbstractVector(0),
              positionsVectorCollection_(),
              typesVector_(0) {

        if ( !particlesVectorVector.empty() ){
            typesVector_ = particlesVectorVector[0].typesVector();
            for (const auto &particlesVector : particlesVectorVector) {
                append(particlesVector);
            }
        }
    }
    ParticlesVectorCollection(const PositionsVectorCollection &positionsVectorCollection,
                              const TypesVector<Type> &typesVector)
            : AbstractVector(positionsVectorCollection.numberOfEntities()),
              positionsVectorCollection_(positionsVectorCollection),
              typesVector_(typesVector) {

        assert(numberOfEntities() == positionsVectorCollection_.numberOfEntities()
               && "The number of entities in ParticlesVectorCollection and PositionsVectorCollection must be equal.");

        assert(positionsVectorCollection.numberOfPositionsEntities() == typesVector_.numberOfEntities()
               && "The number of entities in PositionsVector and TypesVector must be equal.");
    }

    const PositionsVectorCollection &positionsVectorCollection() const {
        return positionsVectorCollection_;
    }

    PositionsVectorCollection &positionsVectorCollection() {
        return positionsVectorCollection_;
    }

    double norm(long i, long j) const{
        return positionsVectorCollection_.norm(i, j);
    }

    ParticlesVector<Type> operator[](long i) const {
        return {positionsVectorCollection_[i],typesVector_};
    }

    const TypesVector<Type>& typesVector() const{
        return typesVector_;
    }

    TypesVector<Type> &typesVector() {
        return typesVector_;
    }

    void insert(const ParticlesVector<Type> &particlesVector, long i) {
        if (typesVector_.numberOfEntities() != 0) {
            assert(typesVector_.asEigenVector() == particlesVector.typesVector().asEigenVector()
            && "Typevectors must be identical.");
        }
        else{
            typesVector_ = particlesVector.typesVector();
        }
        positionsVectorCollection_.insert(particlesVector.positionsVector(), i);
        incrementNumberOfEntities();
    }

    void append(const ParticlesVector<Type> &particlesVector) {
        insert(particlesVector,numberOfEntities());
    }

    void prepend(const ParticlesVector<Type> &particlesVector) {
        insert(particlesVector,0);
    }

    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation) override {
        positionsVectorCollection_.permute(permutation);
        typesVector_.permute(permutation);
    }

    friend std::ostream& operator<<(std::ostream& os, const ParticlesVectorCollection<Type> & pvc){
        for (long i = 0; i < pvc.numberOfEntities(); i++) {

            os << "Vector " <<std::to_string(i + 1) << ":"
               << std::endl
               << pvc[i]
               << std::endl;
        }
        return os;
    }


protected:
    PositionsVectorCollection positionsVectorCollection_;
    TypesVector<Type> typesVector_;
};

using TypedParticlesVectorCollection = ParticlesVectorCollection<int>;
using ElectronsVectorCollection = ParticlesVectorCollection<Spin>;
using AtomsVectorCollection = ParticlesVectorCollection<Element>;


// copy of particlesvector => template?
namespace YAML {
    template<typename Type> struct convert<ParticlesVectorCollection<Type>> {
        static Node encode(const ParticlesVectorCollection<Type> & pv){
            Node node;
            node["Types"] = pv.typesVector();
            node["Positions"] = pv.positionsVectorCollection();
            return node;

        }
        static bool decode(const Node& nodes, ParticlesVectorCollection<Type> & rhs){
            if(!nodes.IsMap())
                return false;
            ParticlesVectorCollection<Type> pv(
                    nodes["Types"].as<TypesVector<Type>>(),
                    nodes["PositionsVectorCollection"].as<PositionsVectorCollection>());
            rhs = pv;
            return true;
        }
    };

    template<typename Type>
    Emitter& operator<< (Emitter& out, const ParticlesVectorCollection<Type>& pv){
        out << BeginMap
            << Key << "Types" << Value << pv.typesVector() << Newline
            << Key << "PositionsVectorCollection" <<  pv.positionsVectorCollection()
            << EndMap;
        return out;
    };
}

#endif //INPSIGHTS_PARTICLESVECTORCOLLECTION_H
