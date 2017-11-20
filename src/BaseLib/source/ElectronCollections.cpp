//
// Created by Michael Heuer on 30.10.17.
//

#include "ElectronCollections.h"

ElectronCollections::ElectronCollections()
        : spinTypeCollection_(SpinTypeCollection())
{}

ElectronCollections::ElectronCollections(const SpinTypeCollection &spinTypeCollection)
        : spinTypeCollection_(spinTypeCollection)
{}



ElectronCollections::ElectronCollections(const std::vector<ElectronCollection> &electronCollections)
        : spinTypeCollection_(static_cast<SpinTypeCollection>(electronCollections[0])) {

    if ( !electronCollections.empty() ){
        for (const auto &electronCollection : electronCollections) {
            this->append(electronCollection);

            // cast the first electronCollection into a SpinTypeCollection
            assert(static_cast<SpinTypeCollection>(electronCollection).spinTypesAsEigenVector() ==
                            spinTypeCollection_.spinTypesAsEigenVector()
                    && "All electron collections must have the same spin type.");
        }
    }
}

ElectronCollections::ElectronCollections(const std::vector<ParticleCollection> &particleCollections,
                                         const SpinTypeCollection &spinTypeCollection)
        : spinTypeCollection_(spinTypeCollection) {


    for (const auto &particleCollection : particleCollections) {
        this->append(ElectronCollection(particleCollection,spinTypeCollection));
    }

    //TODO remove and create ElectronCollections instead => the check is done therein
    /*assert(particleCollections[0].numberOfParticles() == spinTypeCollection_.numberOfSpinTypes()
           && "The number of positions in ParticleCollection and the number of spin type in SpinTypeCollection must match.");

    if ( !particleCollections.empty() ){
        for (const auto &particleCollection : particleCollections) {
            static_cast<ParticleCollections*>(this)->append(particleCollection);

            // cast the first electronCollection into a SpinTypeCollection
            assert(particleCollection.numberOfParticles() == spinTypeCollection_.numberOfSpinTypes()
                   && "The number of positions in ParticleCollection and the number of spin type in SpinTypeCollection must match.");
        }
    }*/
}

ElectronCollections::ElectronCollections(const std::vector<ParticleCollection> &particleCollections)
        : ParticleCollections(particleCollections) {

    if ( !particleCollections.empty() ){
        spinTypeCollection_ = SpinTypeCollection(particleCollections[0].numberOfParticles());
    }
}

void ElectronCollections::insert(const ElectronCollection &electronCollection, long i) {

    // use the SpinTypeCollection if ElectronCollections is empty
    if(getNumberOfParticleCollections() == 0){
        assert( i != 0 && "The collection is empty but the index i is not equal to zero ");
        spinTypeCollection_ = static_cast<SpinTypeCollection>(electronCollection);
    } else {
        assert(spinTypeCollection_.spinTypesAsEigenVector() == electronCollection.spinTypesAsEigenVector()
               && "The spin types of the electron collection to be inserted must match with those already stored.");
    }
    ParticleCollections::insert(electronCollection, i);
}

void ElectronCollections::append(const ElectronCollection &electronCollection) {
    assert(spinTypeCollection_.spinTypesAsEigenVector() == electronCollection.spinTypesAsEigenVector()
           && "The spin types of the electron collection to be inserted must match with those already stored.");
    ParticleCollections::append(electronCollection);
}

void ElectronCollections::prepend(const ElectronCollection &electronCollection) {
    ParticleCollections::insert(electronCollection,0);
}

ElectronCollection ElectronCollections::getElectronCollection(long i) const {
    return ElectronCollection(this->operator[](i),spinTypeCollection_);
}

SpinTypeCollection ElectronCollections::getSpinTypeCollection() const {
    return spinTypeCollection_;
}
