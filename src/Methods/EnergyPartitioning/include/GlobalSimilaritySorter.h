//
// Created by Michael Heuer on 25.09.18.
//

#ifndef AMOLQCPP_GLOBALSIMILARITYSORTER_H
#define AMOLQCPP_GLOBALSIMILARITYSORTER_H

#include "Reference.h"
#include "Sample.h"
#include <Logger.h>
#include <HungarianHelper.h>
#include <spdlog/spdlog.h>
#include <vector>

class GlobalSimilaritySorter {
public:

    GlobalSimilaritySorter(
            std::vector<Reference>& references,
            std::vector<SimilarReferences>& similarReferencesVector,
            double distThresh = 0.1)
            :
            references_(references),
            similarReferencesVector_(similarReferencesVector),
            distThresh_(distThresh),
            console(spdlog::get(Logger::name))
    {
        if(!console){
            Logger::initialize();
            console = spdlog::get(Logger::name);
        };
    }

    bool sort(){
        if(references_.empty()) {
            console->error("References are empty.");
            return false;
        } else if (references_.size() == 1) {
            console->warn("No sorting because only one reference was found.");
            return true;
        }

        auto beginIt = references_.begin();

        if(similarReferencesVector_.empty()) {
            similarReferencesVector_.emplace_back(SimilarReferences(references_.begin()));
            beginIt++;
        }

        for (auto it = beginIt; it != references_.end(); ++it) {
            bool isSimilarQ = false;

            for (auto &simRefs : similarReferencesVector_) {

                auto bestMatch = Metrics::bestMatch<Eigen::Infinity, 2>(
                        (*it).maximum_,
                        (*simRefs.repRefIt_).maximum_);

                if (bestMatch.first < distThresh_) {
                    simRefs.similarReferences_.emplace_back(SimilarReference(it, bestMatch.second));
                    isSimilarQ = true;
                }
            }
            if (!isSimilarQ) similarReferencesVector_.emplace_back(SimilarReferences(it));
        }
        return true;
    }

private:
    std::vector<Reference>& references_;
    std::vector<SimilarReferences>& similarReferencesVector_;
    double distThresh_;
    std::shared_ptr<spdlog::logger> console;
};


#endif //AMOLQCPP_GLOBALSIMILARITYSORTER_H
