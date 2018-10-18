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
            std::vector<Sample>& samples,
            std::vector<Reference>& references,
            std::vector<SimilarReferences>& similarReferencesVector,
            double distThresh = 0.1)
            :
            samples_(samples),
            references_(references),
            similarReferencesVector_(similarReferencesVector),
            distThresh_(distThresh),
            increment_((*references.rbegin()).value() * 1e-4),
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


        for (auto ref = beginIt; ref != references_.end(); ++ref) {
            bool isSimilarQ = false;

            auto simRefsBeginIt = similarReferencesVector_.begin();
            auto simRefsEndIt = similarReferencesVector_.end();
            while (simRefsBeginIt != simRefsEndIt){
                // get bound
                std::vector<Reference> fakeRef = {Reference((*simRefsBeginIt).valueStats().cwiseMax()[0] + increment_)};
                auto simRefUpperBoundIt = std::upper_bound(
                        simRefsBeginIt,similarReferencesVector_.end(),
                        SimilarReferences(fakeRef.begin()));
                //std::cout << std::distance(simRefsBeginIt,simRefUpperBoundIt) << std::endl;

                for (auto simRefs = simRefsBeginIt; simRefs != simRefUpperBoundIt; ++simRefs) {

                    auto bestMatch = Metrics::bestMatch<Eigen::Infinity, 2>(
                            (*ref).maximum(),
                            (*simRefs).representativeReference().maximum());

                    if (bestMatch.first < distThresh_) {
                        (*ref).permute(bestMatch.second, samples_);
                        (*simRefs).add(ref);
                        isSimilarQ = true;
                    }
                }
                if (!isSimilarQ) similarReferencesVector_.emplace_back(SimilarReferences(ref));
                simRefsBeginIt = simRefUpperBoundIt;
            }

        }
        return true;
    }

private:
    std::vector<Sample>& samples_;
    std::vector<Reference>& references_;
    std::vector<SimilarReferences>& similarReferencesVector_;
    double distThresh_, increment_;
    std::shared_ptr<spdlog::logger> console;
};


#endif //AMOLQCPP_GLOBALSIMILARITYSORTER_H
