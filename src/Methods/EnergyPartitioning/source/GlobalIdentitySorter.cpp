//
// Created by heuer on 12.12.18.
//

#include <GlobalIdentitySorter.h>
#include <HungarianHelper.h>

GlobalIdentiySorter::GlobalIdentiySorter(std::vector<Reference> &references, std::vector<Sample> &samples,
                                         double distThresh, double increment)
        :
        references_(references),
        samples_(samples),
        increment_(increment),
        distThresh_(distThresh),
        console(spdlog::get(Logger::name)) {
    if (!console) {
        Logger::initialize();
        console = spdlog::get(Logger::name);
    };
}

bool GlobalIdentiySorter::sort() {
    if (references_.empty()) {
        console->error("References are empty.");
        return false;
    } else if (references_.size() == 1) {
        console->warn("No sorting because only one reference was found.");
        return true;
    }

    std::sort(references_.begin(), references_.end());
    auto beginIt = references_.begin();

    while (beginIt != references_.end()) {
        auto total = std::distance(references_.begin(), references_.end());
        auto endIt = std::upper_bound(beginIt, references_.end(), Reference((*beginIt).value() + increment_));

        console->info("Global identiy search in interval {} to {}, total: {}",
                      total - std::distance(beginIt, references_.end()),
                      total - std::distance(endIt, references_.end()),
                      std::distance(references_.begin(), references_.end()));

        auto it = beginIt;

        if (beginIt != endIt) {
            it++; // start with the element next to beginIt
            while (it != endIt) subLoop(beginIt, it, endIt);
            beginIt = endIt;
        } else ++beginIt; // range is zero
    }
    return true;
}


void GlobalIdentiySorter::subLoop(
        std::vector<Reference>::iterator &beginIt,
        std::vector<Reference>::iterator &it,
        std::vector<Reference>::iterator &endIt) {

    //TODO calculate only alpha electron distances and skip beta electron hungarian if dist is too large
    auto bestMatch = Metrics::spinSpecificBestMatch((*it).maximum(), (*beginIt).maximum());

    if ((*beginIt).maximum().typesVector().multiplicity() == 1) { // consider spin flip

        auto bestMatchFlipped =
                Metrics::spinSpecificBestMatch<Eigen::Infinity, 2>((*it).maximum(), (*beginIt).maximum(), true);

        if ((bestMatch.first <= distThresh_) || (bestMatchFlipped.first <= distThresh_)) {
            if (bestMatch.first <= bestMatchFlipped.first)
                addReference(beginIt, it, bestMatch.second);
            else
                addReference(beginIt, it, bestMatchFlipped.second);
            endIt = std::upper_bound(beginIt, references_.end(), Reference((*beginIt).value() + increment_));
        } else it++;
    } else {  // don't consider spin flip
        if ((bestMatch.first <= distThresh_)) {
            addReference(beginIt, it, bestMatch.second);
            endIt = std::upper_bound(beginIt, references_.end(), Reference((*beginIt).value() + increment_));
        } else it++;
    }
}

// TODO This method should be located inside of a reference container class
void GlobalIdentiySorter::addReference(
        const std::vector<Reference>::iterator &beginIt,
        std::vector<Reference>::iterator &it,
        const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch) const {

    (*it).permute(bestMatch, samples_);
    (*beginIt).mergeReference(it);
    it = references_.erase(it); // erase returns the iterator of the following element
}