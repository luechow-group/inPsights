//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <Sample.h>
#include <Logger.h>
#include <HungarianHelper.h>
#include <algorithm>
#include <DensityBasedScan.h>
#include <CoulombPotential.h>
#include <Statistics.h>

//TODO method header is unclear
double mostDeviatingParticleDistance(
        PositionsVector permutee,
        const Eigen::PermutationMatrix<Eigen::Dynamic> &perm,
        const PositionsVector &ref) {
    assert(permutee.numberOfEntities() == ref.numberOfEntities());

    permutee.permute(perm);
    return Metrics::positionDistancesVector(permutee,ref).lpNorm<Eigen::Infinity>();
}

class GlobalIdentiySorter{
public:

    GlobalIdentiySorter(std::vector<Reference>& references, std::vector<Sample>& samples, double increment, double distThresh)
    :
    references_(references),
    samples_(samples),
    increment_(increment),
    distThresh_(distThresh),
    console(spdlog::get("console"))
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
        auto uit = references_.begin();

        while (lit != references_.end()){
            uit = std::upper_bound(lit,references_.end(),Reference((*lit).negLogSqrdProbabilityDensity_+increment_));
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

    std::vector<Reference>& references_;
    std::vector<Sample>& samples_;
    double increment_,distThresh_;
    std::shared_ptr<spdlog::logger> console;

};


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
    console(spdlog::get("console"))
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
        uit = std::upper_bound(references_.begin(),references_.end(),Reference((*lit).negLogSqrdProbabilityDensity_+increment_));


        while (lit != references_.end()) {
            uit = std::upper_bound(lit,references_.end(),Reference((*lit).negLogSqrdProbabilityDensity_+increment_));

            for (auto it = lit; it != uit;  ++it ) {

                bool isSimilarQ = false;
                for (auto &simRefs : similarReferencesVector_) {
                    // check if refIt is similar to simRefs representative reference

                    // calc perm r->sr so that we can store them in similar reference
                    auto costMatrix = Metrics::positionalDistances(
                            (*it).maximum_.positionsVector(),
                            (*simRefs.representativeReferenceIterator).maximum_.positionsVector());
                    auto bestMatch = Hungarian<double>::findMatching(costMatrix);
                    auto dist = mostDeviatingParticleDistance(
                            (*it).maximum_.positionsVector(),
                            bestMatch,
                            (*simRefs.representativeReferenceIterator).maximum_.positionsVector());

                    //console->info("{}", dist);
                    if (dist < distThresh_) {
                        simRefs.similarReferences_.emplace_back(SimilarReference(it, bestMatch));

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


int main(int argc, char *argv[]) {
    Logger::initialize();
    auto console = spdlog::get(Logger::name);
    double identicalDistThresh = 0.01;
    std::vector<Reference> globallyIdenticalMaxima;
    std::vector<Sample> samples;

    RawDataReader reader(globallyIdenticalMaxima,samples);
    reader.read("raw.bin");
    console->info("number of refs {}",globallyIdenticalMaxima.size());
    auto numberOfElectrons = globallyIdenticalMaxima[0].maximum_.numberOfEntities();


    GlobalIdentiySorter globalIdentiySorter(globallyIdenticalMaxima, samples, identicalDistThresh);
    globalIdentiySorter.sort();
    console->info("finished id sort");
    console->flush();


    console->info("total elems {}",std::distance(globallyIdenticalMaxima.begin(),globallyIdenticalMaxima.end()));


    double similarDistThresh = 0.1;
    std::vector<SimilarReferences> similarReferencesVector;
    GlobalSimilaritySorter globalSimilaritySorter(globallyIdenticalMaxima, similarReferencesVector,similarDistThresh);

    globalSimilaritySorter.sort();

    console->info("total elems {}",similarReferencesVector.size());



    Statistics::RunningStatistics<Eigen::VectorXd> EkinStats;
    Statistics::RunningStatistics<Eigen::MatrixXd> EpotStats;


    Eigen::VectorXd ekin;
    Eigen::MatrixXd epot;
    size_t sampleId, count, simCount, totalCount;

    totalCount = 0;

    for(auto& simRefVector : similarReferencesVector){
        console->info("contained similar refs {}",simRefVector.similarReferences_.size());
        //console->info("representative ref has {} founds", (*i.representativeReferenceIterator).associatedSampleIds_.size());

        simCount = 0;

        EkinStats.reset();
        EpotStats.reset();

        // Representative reference
        count = 1+(*simRefVector.representativeReferenceIterator).associatedSampleIds_.size();
        simCount += count;
        sampleId = (*simRefVector.representativeReferenceIterator).id_;


        ekin = samples[sampleId].kineticEnergies_;
        epot = CoulombPotential::energies(samples[sampleId].sample_);

        console->info("rep ref count {}" ,count);
        EkinStats.add(ekin,unsigned(count));
        EpotStats.add(epot,unsigned(count));

        // Iterate over references being similar to the representative reference.
        for(const auto& simRef : simRefVector.similarReferences_){


            //TODO CONTINUE HERE
            globallyIdenticalMaxima[0];

            // Hier werden irgendwo maxima counts vergessen#
            count = 1+(*simRef.it_).associatedSampleIds_.size();
            simCount += count;
            sampleId = (*simRef.it_).id_;

            auto sampleCopy = samples[sampleId].sample_;
            sampleCopy.permute(simRef.perm_);

            ekin = simRef.perm_ * (samples[sampleId].kineticEnergies_);
            epot = CoulombPotential::energies(sampleCopy);

            EkinStats.add(ekin,unsigned(count));
            EpotStats.add(epot,unsigned(count));
        }
        console->info("sim ref count {}",simCount);
        std::cout <<"mean: ("<<EkinStats.getWeightedSum()<<")\n"<< EkinStats.mean().transpose() << std::endl;
        if(simCount >=2)
            std::cout <<"stdv:\n"<< EkinStats.standardDeviation().transpose() << std::endl<< std::endl;

        totalCount += simCount;
    }
    console->info("overall count {}",totalCount);


    //TODO CHECK PERMUTATION

    // CHECK ADDITIONAL +1

    //TODO DBSCAN

}
