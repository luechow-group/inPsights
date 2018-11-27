//
// Created by Michael Heuer on 26.11.18.
//

#include <AtomsVectorLinkedElectronsVector.h>
#include <Metrics.h>
#include <algorithm>
#include <numeric>

AtomsVectorLinkedElectronsVector::AtomsVectorLinkedElectronsVector(
        std::shared_ptr<AtomsVector> linkedAtomsVector)
        : ElectronsVector(), linkedAtomsVector_(std::move(linkedAtomsVector)) {}

AtomsVectorLinkedElectronsVector::AtomsVectorLinkedElectronsVector(std::shared_ptr<AtomsVector> sharedAtomsVector,
                                                                   ElectronsVector ev)
        : ElectronsVector(std::move(ev)), linkedAtomsVector_(std::move(sharedAtomsVector)) {}

std::vector<long> AtomsVectorLinkedElectronsVector::coreElectronsIndices(long k, double threshold) const {
    std::vector<long> indices;

    for (long i = 0; i < numberOfEntities(); ++i)
        if (Metrics::distance(operator[](i).position(), linkedAtomsVector_->operator[](k).position()) < threshold)
            indices.push_back(i);

    return indices;
}

std::vector<long> AtomsVectorLinkedElectronsVector::coreElectronsIndices(double threshold) const {
    std::vector<long> indices;

    for (long k = 0; k < linkedAtomsVector_->numberOfEntities(); ++k) {
        auto kIndices = coreElectronsIndices(k, threshold);
        std::move(kIndices.begin(), kIndices.end(), std::back_inserter(indices));
    }

    std::sort(indices.begin(), indices.end());
    indices.erase(std::unique( indices.begin(), indices.end() ), indices.end()); // erase duplicates
    return indices;
}

std::vector<long> AtomsVectorLinkedElectronsVector::valenceElectronsIndices(double threshold) const {
    std::vector<long> diff, indices = std::vector<long>(size_t(numberOfEntities()));
    std::iota(indices.begin(), indices.end(), 0);

    auto coreIndices = coreElectronsIndices(threshold);
    std::set_difference(indices.begin(), indices.end(), coreIndices.begin(), coreIndices.end(),
                        std::inserter(diff, diff.begin()));
    return diff;
}
