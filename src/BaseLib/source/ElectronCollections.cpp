//
// Created by Michael Heuer on 30.10.17.
//

#include "ElectronCollections.h"

ElectronCollections::ElectronCollections(const Eigen::VectorXi &spinTypes)
        : spinTypes_(spinTypes) {

    assert(spinTypes_.maxCoeff() <= int(Spin::SpinType::alpha));
    assert(spinTypes_.minCoeff() >= int(Spin::SpinType::beta));
}

ElectronCollections::ElectronCollections(const ElectronCollection &electronCollection)
        : ElectronCollections(std::vector<ElectronCollection>({electronCollection})){}

ElectronCollections::ElectronCollections(const std::vector<ElectronCollection> &electronCollections) {

    for (const auto &electronCollection : electronCollections) {
        this->append(electronCollection);
        assert(electronCollection.getSpinTypes().isApprox(this->spinTypes_)
                && "All electron collections must have the same spin types.");
    }
}

ElectronCollections::ElectronCollections(const std::vector<ParticleCollection> &particleCollections,
                                         const Eigen::VectorXi &spinTypes)
        : ParticleCollections(particleCollections),
          spinTypes_(spinTypes) {

    for (const auto &particleCollection : particleCollections) {
        this->append(ElectronCollection(particleCollection,spinTypes));
    }
}

ElectronCollections::ElectronCollections(const std::vector<ParticleCollection> &particleCollections)
        : ParticleCollections(particleCollections),
          spinTypes_(){

    if ( !particleCollections.empty() ){
        auto numberOfParticles = particleCollections[0].numberOfParticles();
        spinTypes_.resize(numberOfParticles);
        spinTypes_ = Eigen::VectorXi::Zero(numberOfParticles);
    }
}

void ElectronCollections::insert(const ElectronCollection &electronCollection, long i) {

    // if ElectronCollections is empty use the SpinTypeCollection,
    // else compare the spinTypes if the ElectronCollections is not empty
    if(length() == 0){
        assert( i == 0 && "If the collection is empty, the index i must be zero ");
        spinTypes_ = electronCollection.getSpinTypes();
    } else {
        assert(spinTypes_.isApprox(electronCollection.getSpinTypes())
                       && "The spin types of the electron collection to be inserted must match with those already stored.");
    }
    ParticleCollections::insert(static_cast<ParticleCollection>(electronCollection), i);
}

void ElectronCollections::append(const ElectronCollection &electronCollection) {
    this->insert(electronCollection, ParticleCollections::length());
}

void ElectronCollections::prepend(const ElectronCollection &electronCollection) {
    ParticleCollections::insert(electronCollection,0);
}

ElectronCollection ElectronCollections::getElectronCollection(long i) const {
    return ElectronCollection(this->operator[](i),spinTypes_);
}

Eigen::VectorXi ElectronCollections::getSpinTypes() const {
    return spinTypes_;
}
