//
// Created by Michael Heuer on 30.10.17.
//

#ifndef AMOLQCGUI_ELECTRONCOLLECTIONPATH_H
#define AMOLQCGUI_ELECTRONCOLLECTIONPATH_H

#include "ParticlesVectorCollection.h"
#include "ElectronsVector.h"
#include "SpinTypesVector.h"

class ElectronsVectorCollection : public ParticlesVectorCollection{
public:
    ElectronsVectorCollection();
    explicit ElectronsVectorCollection(const SpinTypesVector& spinTypesVector);
    explicit ElectronsVectorCollection(const ElectronsVector& electronsVector);
    explicit ElectronsVectorCollection(const std::vector<ElectronsVector>& electronsVectorVector);
    explicit ElectronsVectorCollection(const PositionsVectorCollection& electronsVector);

    explicit ElectronsVectorCollection(const PositionsVectorCollection& electronsVector,
                                 const SpinTypesVector& spinTypesVector);

    ElectronsVector operator[](long i) const;

    const SpinTypesVector& spinTypesVector() const;
    SpinTypesVector& spinTypesVector();

    void insert (const ElectronsVector& electronsVector, long i);
    void append (const ElectronsVector& electronsVector);
    void prepend(const ElectronsVector& electronsVector);
    void permute(long i, long j) override;

private:
    SpinTypesVector spinTypesVector_;
};

#endif //AMOLQCGUI_ELECTRONCOLLECTIONPATH_H
