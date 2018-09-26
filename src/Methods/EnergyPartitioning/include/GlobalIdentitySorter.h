//
// Created by Michael Heuer on 25.09.18.
//

#ifndef AMOLQCPP_GLOBALIDENTITYSORTER_H
#define AMOLQCPP_GLOBALIDENTITYSORTER_H


#include "Reference.h"
#include "Sample.h"
#include <Logger.h>
#include <HungarianHelper.h>
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
    {}

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
        auto lit = references_.begin();

        while (lit != references_.end()){
            auto uit = std::upper_bound(lit,references_.end(),Reference((*lit).negLogSqrdProbabilityDensity_+increment_));
            auto it = lit;
            it++; // start with the element next to lit

            while(it != uit) subLoop(lit,it,uit);
            lit = uit;
        }
        return true;
    }

private:
    void subLoop(
            std::vector<Reference>::iterator& lit,
            std::vector<Reference>::iterator& it,
            std::vector<Reference>::iterator& uit) {

        //TODO calculate only alpha electron distances and skip beta electron hungarian if dist is too large
        auto bestMatch = Metrics::spinSpecificBestMatchNorm((*it).maximum_, (*lit).maximum_);

        if((*lit).maximum_.typesVector().multiplicity() == 1) { // consider spin flip
            
            auto bestMatchFlipped = Metrics::spinSpecificBestMatchNorm<Eigen::Infinity,2>((*it).maximum_, (*lit).maximum_, true);

            if( (bestMatch.first <= distThresh_) || (bestMatchFlipped.first <= distThresh_) ){
                if(bestMatch.first <= bestMatchFlipped.first)
                    addReference(lit, it, bestMatch.second);
                else
                    addReference(lit, it, bestMatchFlipped.second);
                uit = std::upper_bound(lit, references_.end(), Reference((*lit).negLogSqrdProbabilityDensity_+increment_));
            }
            else it++;
        }
        else {  // don't consider spin flip
            if( (bestMatch.first <= distThresh_) ) {
                addReference(lit, it, bestMatch.second);
                uit = std::upper_bound(lit,references_.end(),Reference((*lit).negLogSqrdProbabilityDensity_+increment_));
            }
            else it++;
        }
    }

    //TODO find better name
    void addReference(const std::vector<Reference, std::allocator<Reference>>::iterator &lit,
            std::vector<Reference, std::allocator<Reference>>::iterator &it,
    const Eigen::PermutationMatrix<Eigen::Dynamic> &bestMatch) const {
        samples_[(*it).id_].sample_.permute(bestMatch);

        (*lit).addAssociation((*it).id_);
        (*lit).associatedSampleIds_.insert(
                (*lit).associatedSampleIds_.end(),
                make_move_iterator((*it).associatedSampleIds_.begin()),
                make_move_iterator((*it).associatedSampleIds_.end())
        );

        it = references_.erase(it); // returns the iterator of the following element
    }

    std::vector<Reference>& references_;
    std::vector<Sample>& samples_;
    double increment_,distThresh_;
    std::shared_ptr<spdlog::logger> console;

};


#endif //AMOLQCPP_GLOBALIDENTITYSORTER_H
