//
// Created by Michael Heuer on 26.11.18.
//

#ifndef INPSIGHTS_ATOMSVECTORLINKEDELECTRONSVECTOR_H
#define INPSIGHTS_ATOMSVECTORLINKEDELECTRONSVECTOR_H

#include <ParticlesVector.h>
#include <memory>

class AtomsVectorLinkedElectronsVector : public ElectronsVector {
public:

    explicit AtomsVectorLinkedElectronsVector(std::shared_ptr<AtomsVector> linkedAtomsVector);

    AtomsVectorLinkedElectronsVector(std::shared_ptr<AtomsVector> sharedAtomsVector, ElectronsVector ev);

    std::vector<long> coreElectronsIndices(long k, double threshold = 0.1) const;

    std::vector<long> coreElectronsIndices(double threshold = 0.1) const;

    std::vector<long> valenceElectronsIndices(double threshold = 0.1) const;

private:
    std::shared_ptr<AtomsVector> linkedAtomsVector_;
};

#endif //INPSIGHTS_ATOMSVECTORLINKEDELECTRONSVECTOR_H
