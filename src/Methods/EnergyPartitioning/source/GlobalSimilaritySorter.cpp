//
// Created by heuer on 12.12.18.
//

#include <GlobalSimilaritySorter.h>
#include <ValueSorter.h>

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
// assumes a sorted reference vector
bool GlobalSimilaritySorter::sort() {
    // first, sort references by value
    ValueSorter::sortReferencesByValue(references_);

    // insert first element
    similarReferencesVector_.emplace_back(SimilarReferences(references_.begin()));
    // start with the second reference
    for (auto ref = references_.begin()+1; ref != references_.end(); ++ref) {
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
    }

    return true;
}