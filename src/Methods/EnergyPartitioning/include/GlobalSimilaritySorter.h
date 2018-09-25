//
// Created by Michael Heuer on 25.09.18.
//

#ifndef AMOLQCPP_GLOBALSIMILARITYSORTER_H
#define AMOLQCPP_GLOBALSIMILARITYSORTER_H


class GlobalSimilaritySorter {
public:

    GlobalSimilaritySorter(
            std::vector<Reference>& references,
            std::vector<SimilarReferences>& similarReferencesVector,
            double increment, double distThresh)
            :
            references_(references),
            similarReferencesVector_(similarReferencesVector),
            increment_(increment),
            distThresh_(distThresh),
            console(spdlog::get(Logger::name))
    {}

    GlobalSimilaritySorter(
            std::vector<Reference>& references,
            std::vector<SimilarReferences>& similarReferencesVector,
            double distThresh = 0.1)
            :
            GlobalSimilaritySorter(references, similarReferencesVector,
                                   (*references.rbegin()).negLogSqrdProbabilityDensity_ * 1e-5, distThresh)
    {}


    bool sort(){

        if(similarReferencesVector_.empty())
            similarReferencesVector_.emplace_back(SimilarReferences(references_.begin()));

        // for all references in range
        // for (auto refIt = references_.begin()+1; refIt != references_.end(); ++refIt) {
        auto lit = references_.begin();
        auto uit = references_.begin();

        while (lit != references_.end()) {
            uit = std::upper_bound(lit,references_.end(), Reference((*lit).negLogSqrdProbabilityDensity_+increment_));

            for (auto it = lit; it != uit;  ++it ) {

                bool isSimilarQ = false;
                for (auto &simRefs : similarReferencesVector_) {
                    // check if refIt is similar to simRefs representative reference

                    // calc perm r->sr so that we can store them in similar reference
                    auto bestMatch = Metrics::bestMatchNorm<Eigen::Infinity>(
                            (*it).maximum_.positionsVector(),
                            (*simRefs.representativeReferenceIterator).maximum_.positionsVector());

                    console->info("{}", bestMatch.first);
                    if (bestMatch.first < distThresh_) {
                        simRefs.similarReferences_.emplace_back(SimilarReference(it, bestMatch.second));

                        console->info("merged");
                        isSimilarQ = true;
                    }
                }
                if (!isSimilarQ) similarReferencesVector_.emplace_back(SimilarReferences(it));
            }
            lit = uit;
        }
        return true;
    }
private:

    std::vector<Reference>& references_;
    std::vector<SimilarReferences>& similarReferencesVector_;
    double increment_,distThresh_;
    std::shared_ptr<spdlog::logger> console;
};


#endif //AMOLQCPP_GLOBALSIMILARITYSORTER_H
