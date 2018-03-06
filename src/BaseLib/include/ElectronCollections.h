//
// Created by Michael Heuer on 30.10.17.
//

#ifndef AMOLQCGUI_ELECTRONCOLLECTIONPATH_H
#define AMOLQCGUI_ELECTRONCOLLECTIONPATH_H

#include "ParticleCollections.h"
#include "ElectronCollection.h"

class ElectronCollections : public ParticleCollections{
public:
    explicit ElectronCollections(const Eigen::VectorXi &spinTypes);
    explicit ElectronCollections(const ElectronCollection &electronCollection);
    explicit ElectronCollections(const std::vector<ElectronCollection> &electronCollections);
    explicit ElectronCollections(const std::vector<ParticleCollection> &particleCollections);
    explicit ElectronCollections(const std::vector<ParticleCollection> &particleCollections,
                                 const Eigen::VectorXi &spinTypes);

    ElectronCollection getElectronCollection(long i) const;

    void insert (const ElectronCollection& electronCollection, long i);
    virtual void append (const ElectronCollection& electronCollection);
    void prepend(const ElectronCollection& electronCollection);

    Eigen::VectorXi getSpinTypes() const;

private:
    Eigen::VectorXi spinTypes_;
};

#endif //AMOLQCGUI_ELECTRONCOLLECTIONPATH_H
