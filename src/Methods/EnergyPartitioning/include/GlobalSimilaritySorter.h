//
// Created by Michael Heuer on 25.09.18.
//

#ifndef INPSIGHTS_GLOBALSIMILARITYSORTER_H
#define INPSIGHTS_GLOBALSIMILARITYSORTER_H

#include "SimilarReferences.h"
#include <Logger.h>
#include <HungarianHelper.h>
#include <spdlog/spdlog.h>
#include <vector>

class GlobalSimilaritySorter {
public:

    GlobalSimilaritySorter(std::vector<Sample> &samples, std::vector<Reference> &references,
                               std::vector<SimilarReferences> &similarReferencesVector,
                               double distThresh, double increment);

    GlobalSimilaritySorter(
            std::vector<Sample>& samples,
            std::vector<Reference>& references,
            std::vector<SimilarReferences>& similarReferencesVector,
            double distThresh = 0.1);

    bool sort();

private:
    std::vector<Sample>& samples_;
    std::vector<Reference>& references_;
    std::vector<SimilarReferences>& similarReferencesVector_;
    double distThresh_, increment_;
    std::shared_ptr<spdlog::logger> console;
};


#endif //INPSIGHTS_GLOBALSIMILARITYSORTER_H
