//
// Created by heuer on 12.12.18.
//

#include <GlobalSimilaritySorter.h>

GlobalSimilaritySorter::GlobalSimilaritySorter(std::vector<Sample> &samples, std::vector<Reference> &references,
                                               std::vector<SimilarReferences> &similarReferencesVector,
                                               double distThresh, double increment)
        :
        samples_(samples),
        references_(references),
        similarReferencesVector_(similarReferencesVector),
        distThresh_(distThresh),
        increment_(increment),
        console(spdlog::get(Logger::name)) {
    if (!console) {
        Logger::initialize();
        console = spdlog::get(Logger::name);
    };
}

bool GlobalSimilaritySorter::sort() {
    if (references_.empty()) {
        console->error("References are empty.");
        return false;
    } else if (references_.size() == 1) {
        console->warn("No sorting because only one reference was found.");
        return true;
    }

    auto beginIt = references_.begin();

    if (similarReferencesVector_.empty()) {
        similarReferencesVector_.emplace_back(SimilarReferences(references_.begin()));
        beginIt++;
    }

    for (auto ref = beginIt; ref != references_.end(); ++ref) {
        bool isSimilarQ = false;

        std::vector<Reference> lowerRef = {Reference((*ref).value() - increment_)};
        std::vector<Reference> upperRef = {Reference((*ref).value() + increment_)};

        auto simRefLowerBoundIt = std::lower_bound(
                similarReferencesVector_.begin(),
                similarReferencesVector_.end(),
                SimilarReferences(lowerRef.begin()));
        auto simRefUpperBoundIt = std::upper_bound(
                similarReferencesVector_.begin(),
                similarReferencesVector_.end(),
                SimilarReferences(upperRef.begin()));

        for (auto simRefs = simRefLowerBoundIt; simRefs != simRefUpperBoundIt; ++simRefs) {
            auto bestMatch = Metrics::bestMatch<Eigen::Infinity, 2>(
                    (*ref).maximum(),
                    (*simRefs).representativeReference().maximum());

            if (bestMatch.first < distThresh_) {
                (*ref).permute(bestMatch.second, samples_);
                (*simRefs).add(ref);
                isSimilarQ = true;
                break;
            }
        }
        if (!isSimilarQ) {
            similarReferencesVector_.emplace_back(SimilarReferences(ref));
        }
        std::sort(similarReferencesVector_.begin(), similarReferencesVector_.end());
    }
    return true;
}