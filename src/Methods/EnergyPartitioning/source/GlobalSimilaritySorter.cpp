//
// Created by heuer on 12.12.18.
//

#include <GlobalIdentitySorter.h>
#include <GlobalSimilaritySorter.h>
#include <ValueSorter.h>

namespace Settings {
    GlobalSimilaritySorter::GlobalSimilaritySorter() {
        similarityRadius.onChange().connect([&](double value) {
            if(value < ::GlobalIdentitySorter::settings.identityRadius.get())
                throw std::invalid_argument(
                        "The " + similarityRadius.name() + " with " + std::to_string(similarityRadius.get())
                        + " is smaller than the "+ ::GlobalIdentitySorter::settings.identityRadius.name() + " with "
                        + std::to_string(::GlobalIdentitySorter::settings.identityRadius.get()));
        });
    }

    GlobalSimilaritySorter::GlobalSimilaritySorter(const YAML::Node &node)
            : GlobalSimilaritySorter() {
        doubleProperty::decode(node[className], similarityRadius);
        doubleProperty::decode(node[className], valueIncrement);
    }

    void GlobalSimilaritySorter::addToNode(YAML::Node &node) const {
        node[className][similarityRadius.name()] = similarityRadius.get();
        node[className][valueIncrement.name()] = valueIncrement.get();
    }
}
YAML_SETTINGS_DEFINITION(Settings::GlobalSimilaritySorter)

GlobalSimilaritySorter::GlobalSimilaritySorter(
        std::vector<Sample> &samples,
        std::vector<Reference> &references,
        std::vector<SimilarReferences> &similarReferencesVector)
        :
        samples_(samples),
        references_(references),
        similarReferencesVector_(similarReferencesVector) {}
        
// assumes a sorted reference vector
bool GlobalSimilaritySorter::sort() {
    auto similarityRadius = settings.similarityRadius.get();
    auto valueIncrement = settings.valueIncrement.get();

    // first, sort references by value
    ValueSorter::sortReferencesByValue(references_);

    // insert first element
    similarReferencesVector_.emplace_back(SimilarReferences(references_.begin()));
    // start with the second reference
    for (auto ref = references_.begin()+1; ref != references_.end(); ++ref) {
        bool isSimilarQ = false;

        std::vector<Reference> lowerRef = {Reference((*ref).value() - valueIncrement)};
        std::vector<Reference> upperRef = {Reference((*ref).value() + valueIncrement)};

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

            if (bestMatch.first < similarityRadius) {
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