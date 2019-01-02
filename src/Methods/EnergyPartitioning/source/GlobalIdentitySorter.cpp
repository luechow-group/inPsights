//
// Created by heuer on 12.12.18.
//

#include <GlobalIdentitySorter.h>
#include <GlobalSimilaritySorter.h>
#include <HungarianHelper.h>
#include <ValueSorter.h>

namespace Settings {
    GlobalIdentitySorter::GlobalIdentitySorter() {
        identityRadius.onChange().connect([&](double value) {
            if(value > ::GlobalSimilaritySorter::settings.similarityRadius.get())
                throw std::invalid_argument(
                        "The " + identityRadius.name() + " with " + std::to_string(identityRadius.get())
                        + " is greater than the "+ ::GlobalSimilaritySorter::settings.similarityRadius.name() + " with "
                        + std::to_string(::GlobalSimilaritySorter::settings.similarityRadius.get()));
        });
    }

    GlobalIdentitySorter::GlobalIdentitySorter(const YAML::Node &node)
            : GlobalIdentitySorter() {
        doubleProperty::decode(node[className], identityRadius);
        doubleProperty::decode(node[className], identityValueIncrement);
    }

    void GlobalIdentitySorter::appendToNode(YAML::Node &node) const {
        node[className][identityRadius.name()] = identityRadius.get();
        node[className][identityValueIncrement.name()] = identityValueIncrement.get();
    }
}
YAML_SETTINGS_DEFINITION(Settings::GlobalIdentitySorter)

GlobalIdentitySorter::GlobalIdentitySorter(
        std::vector<Reference> &references,
        std::vector<Sample> &samples)
        :
        references_(references),
        samples_(samples) {}


bool GlobalIdentitySorter::sort() {
    auto identityRadius = settings.identityRadius.get();
    auto valueIncrement = settings.identityValueIncrement.get();

    // first, sort references by value
    ValueSorter::sortReferencesByValue(references_);

    auto beginIt = references_.begin();

    while (beginIt != references_.end()) {
        auto total = std::distance(references_.begin(), references_.end());
        auto endIt = std::upper_bound(beginIt, references_.end(), Reference((*beginIt).value() + valueIncrement));

        Logger::console->info("Global identiy search in interval {} to {}, total: {}",
                      total - std::distance(beginIt, references_.end()),
                      total - std::distance(endIt, references_.end()),
                      std::distance(references_.begin(), references_.end()));

        auto it = beginIt;

        if (beginIt != endIt) {
            it++; // start with the element next to beginIt
            while (it != endIt)
                subLoop(beginIt, it, endIt, identityRadius, valueIncrement);

            beginIt = endIt;
        } else ++beginIt; // range is zero
    }
    return true;
}


void GlobalIdentitySorter::subLoop(
        std::vector<Reference>::iterator &beginIt,
        std::vector<Reference>::iterator &it,
        std::vector<Reference>::iterator &endIt,
        double distThresh,
        double valueIncrement) {

    //TODO calculate only alpha electron distances and skip beta electron hungarian if dist is too large
    auto bestMatch = Metrics::spinSpecificBestMatch((*it).maximum(), (*beginIt).maximum());

    if ((*beginIt).maximum().typesVector().multiplicity() == 1) { // consider spin flip

        auto bestMatchFlipped =
                Metrics::spinSpecificBestMatch<Eigen::Infinity, 2>((*it).maximum(), (*beginIt).maximum(), true);

        if ((bestMatch.first <= distThresh) || (bestMatchFlipped.first <= distThresh)) {
            if (bestMatch.first <= bestMatchFlipped.first)
                addReference(beginIt, it, bestMatch.second);
            else
                addReference(beginIt, it, bestMatchFlipped.second);
            endIt = std::upper_bound(beginIt, references_.end(), Reference((*beginIt).value() + valueIncrement));
        } else it++;
    } else {  // don't consider spin flip
        if ((bestMatch.first <= distThresh)) {
            addReference(beginIt, it, bestMatch.second);
            endIt = std::upper_bound(beginIt, references_.end(), Reference((*beginIt).value() + valueIncrement));
        } else it++;
    }
}

// TODO This method should be located inside of a reference container class
void GlobalIdentitySorter::addReference(
        const std::vector<Reference>::iterator &beginIt,
        std::vector<Reference>::iterator &it,
        const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch) const {

    (*it).permute(bestMatch, samples_);
    (*beginIt).mergeReference(it);
    it = references_.erase(it); // erase returns the iterator of the following element
}