//
// Created by Michael Heuer on 25.09.18.
//

#ifndef AMOLQCPP_GLOBALIDENTITYSORTER_H
#define AMOLQCPP_GLOBALIDENTITYSORTER_H


#include "Reference.h"
#include "Sample.h"
#include <Logger.h>
#include <HungarianHelper.h>
#include <spdlog/spdlog.h>
#include <vector>

class GlobalIdentiySorter{
public:

    GlobalIdentiySorter(std::vector<Reference>& references, std::vector<Sample>& samples, double increment, double distThresh)
            :
            references_(references),
            samples_(samples),
            increment_(increment),
            distThresh_(distThresh),
            console(spdlog::get(Logger::name))
    {
        if(!console){
            Logger::initialize();
            console = spdlog::get(Logger::name);
        };
    }

    GlobalIdentiySorter(std::vector<Reference>& references, std::vector<Sample>& samples, double distThresh = 0.01)
            :
            GlobalIdentiySorter(references, samples, (*references.rbegin()).negLogSqrdProbabilityDensity_ * 1e-7, distThresh)
    {}

    bool sort(){
        if(references_.empty()) {
            console->error("References are empty.");
            return false;
        }
        else if (references_.size() == 1) {
            console->warn("No sorting because only one reference was found.");
            return true; // no sorting
        }
        std::sort(references_.begin(),references_.end());
        auto beginIt = references_.begin();

        while (beginIt != references_.end()){
            auto endIt = std::upper_bound(beginIt,references_.end(),Reference((*beginIt).negLogSqrdProbabilityDensity_+increment_));
            auto it = beginIt;
            it++; // start with the element next to beginIt

            while(it != endIt) subLoop(beginIt,it,endIt);
            beginIt = endIt;
        }
        return true;
    }

private:
    void subLoop(
            std::vector<Reference>::iterator& beginIt,
            std::vector<Reference>::iterator& it,
            std::vector<Reference>::iterator& endIt) {

        //TODO calculate only alpha electron distances and skip beta electron hungarian if dist is too large
        auto bestMatch = Metrics::spinSpecificBestMatchNorm((*it).maximum_, (*beginIt).maximum_);

        if((*beginIt).maximum_.typesVector().multiplicity() == 1) { // consider spin flip

            auto bestMatchFlipped = Metrics::spinSpecificBestMatchNorm<Eigen::Infinity,2>((*it).maximum_, (*beginIt).maximum_, true);

            if( (bestMatch.first <= distThresh_) || (bestMatchFlipped.first <= distThresh_) ){
                if(bestMatch.first <= bestMatchFlipped.first)
                    addReference(beginIt, it, bestMatch.second);
                else
                    addReference(beginIt, it, bestMatchFlipped.second);
                endIt = std::upper_bound(beginIt, references_.end(), Reference((*beginIt).negLogSqrdProbabilityDensity_+increment_));
            }
            else it++;
        }
        else {  // don't consider spin flip
            if( (bestMatch.first <= distThresh_) ) {
                addReference(beginIt, it, bestMatch.second);
                endIt = std::upper_bound(beginIt,references_.end(),Reference((*beginIt).negLogSqrdProbabilityDensity_+increment_));
            }
            else it++;
        }
    }

    void addReference(
            const std::vector<Reference, std::allocator<Reference>>::iterator &beginIt,
            std::vector<Reference, std::allocator<Reference>>::iterator &it,
            const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch) const {

        samples_[(*it).id_].sample_.permute(bestMatch);

        (*beginIt).addAssociation((*it).id_);
        (*beginIt).associatedSampleIds_.insert(
                (*beginIt).associatedSampleIds_.end(),
                make_move_iterator((*it).associatedSampleIds_.begin()),
                make_move_iterator((*it).associatedSampleIds_.end())
        );

        it = references_.erase(it); // erase returns the iterator of the following element
    }

    std::vector<Reference>& references_;
    std::vector<Sample>& samples_;
    double increment_,distThresh_;
    std::shared_ptr<spdlog::logger> console;

};


#endif //AMOLQCPP_GLOBALIDENTITYSORTER_H
