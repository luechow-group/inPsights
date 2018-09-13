//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <Sample.h>
#include <Logger.h>
#include <HungarianHelper.h>
#include <algorithm>
#include <DensityBasedScan.h>

//TODO method header is unclear
double mostDeviatingParticleDistance(
        PositionsVector permutee,
        const Eigen::PermutationMatrix<Eigen::Dynamic> &perm,
        const PositionsVector &ref) {
    assert(permutee.numberOfEntities() == ref.numberOfEntities());

    permutee.permute(perm);
    return Metrics::positionDistancesVector(permutee,ref).lpNorm<Eigen::Infinity>();
}

class AbstractSorter{
public:
    AbstractSorter(std::vector<Reference>& references, std::vector<Sample>& samples)
    :
    references_(references),
    samples_(samples),
    console(spdlog::get("console"))
    {}

protected:
    std::vector<Reference>& references_;
    std::vector<Sample>& samples_;
    std::shared_ptr<spdlog::logger> console;
};

class GlobalIdentiySorter : public AbstractSorter{
public:

    GlobalIdentiySorter(std::vector<Reference>& references, std::vector<Sample>& samples, double increment, double distThresh)
    :
    AbstractSorter(references,samples),
    increment_(increment),
    distThresh_(distThresh)
    {}

    GlobalIdentiySorter(std::vector<Reference>& references, std::vector<Sample>& samples, double distThresh = 0.01)
    :
    GlobalIdentiySorter(references, samples, (*references.rbegin()).negLogSqrdProbabilityDensity_ * 1e-7, distThresh)
    {}


    bool doSort(){
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
        auto uit = references_.begin();

        while (lit != references_.end()){
            uit = std::upper_bound(references_.begin(),references_.end(),Reference((*lit).negLogSqrdProbabilityDensity_+increment_));
            //uit = references_.upper_bound(Reference((*lit).negLogSqrdProbabilityDensity_+increment_));
            auto it = lit;
            it++; // start with the element next to lit

            while(it != uit) subLoop(lit,it);

            lit = uit;
        }

        return true;
    }

private:
    void subLoop(std::vector<Reference>::iterator& lit, std::vector<Reference>::iterator& it){

        /*PUT INTO METHOD*/
        //TODO CHECK MULTIPLICITY
        //TODO maybe calculate only alpha electron distance and skip beta electron hungarian if dist is too high already
        auto bestMatch = HungarianHelper::spinSpecificHungarian((*it).maximum_,(*lit).maximum_);
        auto bestMatchFlip = HungarianHelper::spinSpecificHungarian((*it).maximum_,(*lit).maximum_,true);

        double dist= mostDeviatingParticleDistance(
                (*it).maximum_.positionsVector(), bestMatch,
                (*lit).maximum_.positionsVector());

        double distFlip = mostDeviatingParticleDistance(
                (*it).maximum_.positionsVector(), bestMatchFlip,
                (*lit).maximum_.positionsVector());
        /*PUT INTO METHOD END*/

        console->info("{},{}",dist,distFlip);
        if( (dist <= distThresh_) || (distFlip <= distThresh_) ){
            // refs are identical

            console->info("MERGE");

            if(dist <= distFlip)
                samples_[(*it).id_].sample_.permute(bestMatch);
            else
                samples_[(*it).id_].sample_.permute(bestMatchFlip);

            (*lit).addAssociation((*it).id_);
            //(*lit).associatedSamples_.merge((*it).associatedSamples_);
            (*lit).associatedSampleIds_.insert(
                    (*lit).associatedSampleIds_.end(),
                    std::make_move_iterator((*it).associatedSampleIds_.begin()),
                    std::make_move_iterator((*it).associatedSampleIds_.end())
                    );

            it = references_.erase(it); // returns the iterator of the following element
        } else {
            it++;
        }
    }

    double increment_,distThresh_;
};

int main(int argc, char *argv[]) {
    Logger::initialize();
    auto console = spdlog::get(Logger::name);
    double identicalDistThresh = 0.01;

    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;

    RawDataReader reader(globallyIdenticalMaxima,samples);
    reader.read("raw.bin");

    console->info("number of refs {}",globallyIdenticalMaxima.size());

    GlobalIdentiySorter globalIdentiySorter(globallyIdenticalMaxima, samples, identicalDistThresh);
    globalIdentiySorter.doSort();

    console->info("finished id sort");
    console->flush();


    std::vector<SimilarReferencesCollection> similarReferencesCollections;
    // add first ref

    console->info("total elems {}",std::distance(globallyIdenticalMaxima.begin(),globallyIdenticalMaxima.end()));





    //TODO add lower + upper bound
    for (auto it = globallyIdenticalMaxima.begin(); it != globallyIdenticalMaxima.end();  ++it ){
        console->info("elem {}",std::distance(globallyIdenticalMaxima.begin(),it));
        // check for similarity
        bool isSimilarQ = false;
        for(auto& sr : similarReferencesCollections){

            // calc perm r->sr so that we can store them in sr
            auto costMatrix = Metrics::positionalDistances(
                    (*it).maximum_.positionsVector(),
                    (*sr.representativeReferenceIterator).maximum_.positionsVector());
            auto bestMatch = Hungarian<double>::findMatching(costMatrix);
            auto dist = mostDeviatingParticleDistance(
                    (*it).maximum_.positionsVector(),
                    bestMatch,
                    (*sr.representativeReferenceIterator).maximum_.positionsVector());

            console->info("{}",dist);
            if(dist < identicalDistThresh) {
                sr.similarReferences_.emplace_back(SimilarReference(it, bestMatch));

                isSimilarQ = true;
            }
        }
        if(!isSimilarQ) similarReferencesCollections.emplace_back(SimilarReferencesCollection(it));
    }


    //TODO add DBSCAN

}
