//
// Created by heuer on 11.04.19.
//

#include <Sample.h>

Sample::Sample(ElectronsVector sample, Eigen::VectorXd kineticEnergies)
: sample_(std::move(sample)), kineticEnergies_(std::move(kineticEnergies)) {
    assert(sample.numberOfEntities() == kineticEnergies_.size());
}

void Sample::permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& perm) {
    assert(perm.indices().size() == sample_.numberOfEntities()
           && "Permutation vector size must match with the number of entities");
    sample_.permute(perm);
    kineticEnergies_ = perm*kineticEnergies_;
}
