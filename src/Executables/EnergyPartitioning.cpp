//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <Sample.h>
#include <Logger.h>
#include <HungarianHelper.h>

double mostDeviatingParticleDistance(const PositionsVector& ref, PositionsVector other,
        const Eigen::PermutationMatrix<Eigen::Dynamic>& perm) {
    assert(ref.numberOfEntities() == other.numberOfEntities());

    other.permute(perm);
    return Metrics::positionDistancesVector(ref,other).lpNorm<Eigen::Infinity>();
}

class AbstractSorter{
public:
    AbstractSorter(std::set<Reference>& references, std::vector<Sample>& samples)
    :
    references_(references),
    samples_(samples),
    console(spdlog::get("console"))
    {}

protected:
    std::set<Reference>& references_;
    std::vector<Sample>& samples_;
    std::shared_ptr<spdlog::logger> console;
};

class GlobalIdentiySorter : public AbstractSorter{
public:

    GlobalIdentiySorter(std::set<Reference>& references, std::vector<Sample>& samples, double increment, double distThresh)
    :
    AbstractSorter(references,samples),
    increment_(increment),
    distThresh_(distThresh)
    {}

    GlobalIdentiySorter(std::set<Reference>& references, std::vector<Sample>& samples, double distThresh = 0.01)
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

        auto lit = references_.begin();
        auto uit = references_.begin();

        while (lit != references_.end()){
            uit = references_.upper_bound(Reference((*lit).negLogSqrdProbabilityDensity_+increment_));
            auto it = lit;
            it++; // start with the element next to lit

            while(it != uit) subLoop(lit,it);

            lit = uit;
        }

        return true;
    }

private:
    void subLoop(std::set<Reference>::iterator& lit, std::set<Reference>::iterator& it){

        //TODO CHECK MULTIPLICITY
        auto bestMatch = HungarianHelper::spinSpecificHungarian((*it).maximum_,(*lit).maximum_);
        auto bestMatchFlip = HungarianHelper::spinSpecificHungarian((*it).maximum_,(*lit).maximum_,true);

        double dist= mostDeviatingParticleDistance(
                (*lit).maximum_.positionsVector(),
                (*it).maximum_.positionsVector(),bestMatch);

        double distFlip = mostDeviatingParticleDistance(
                (*lit).maximum_.positionsVector(),
                (*it).maximum_.positionsVector(),bestMatchFlip);

        if( (dist <= distThresh_) || (distFlip <= distThresh_) ){
            // refs are identical

            if(dist <= distFlip)
                samples_[(*it).id_].sample_.permute(bestMatch);
            else
                samples_[(*it).id_].sample_.permute(bestMatchFlip);

            (*lit).addAssociation((*it).id_);
            (*lit).associations_.merge((*it).associations_);

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

    std::set<Reference> references;
    std::vector<Sample> samples;

    RawDataReader reader(references,samples);
    reader.read("raw.bin");

    console->info("number of refs {}",references.size());

    GlobalIdentiySorter globalIdentiySorter(references,samples);
    globalIdentiySorter.doSort();

    console->flush();

    for (auto a : references){
        std::cout << a.id_ << " {";
        for (auto b : a.associations_){
            std::cout << b << " ";
        }
        std::cout << " }" << std::endl;
    }

}
/*console->info("  it={}", (std::distance(references.begin(),it)) );
            // check hungarian

            //TODO CHECK MULTIPLICITY

            auto bestMatch = HungarianHelper::spinSpecificHungarian((*it).maximum_,(*lit).maximum_);
            auto bestMatchFlip = HungarianHelper::spinSpecificHungarian((*it).maximum_,(*lit).maximum_,true);

            double dist= mostDeviatingParticleDistance(
                    (*lit).maximum_.positionsVector(),
                    (*it).maximum_.positionsVector(),bestMatch);

            double distFlip = mostDeviatingParticleDistance(
                    (*lit).maximum_.positionsVector(),
                    (*it).maximum_.positionsVector(),bestMatchFlip);

            console->info("dist: {} {}",dist,distFlip);

            if( (dist <= distThresh) || (distFlip <= distThresh) ){
                // refs are identical

                if(dist <= distFlip)
                    samples[(*it).id_].sample_.permute(bestMatch);
                else
                    samples[(*it).id_].sample_.permute(bestMatchFlip);

                console->info("uit={}", (std::distance(references.begin(),uit)) );
                console->info("MERGE {}<-{}",(*lit).id_,(*it).id_);


                (*lit).addAssociation((*it).id_);
                (*lit).associations_.merge((*it).associations_);

                it = references.erase(it); // returns the iterator of the following element
                console->info("uit={}", (std::distance(references.begin(),uit)) );

            } else {
                it++;
            }

            console->info("end  it={}", (std::distance(references.begin(),it)) );
            console->info("end uit={}", (std::distance(references.begin(),uit)) );*/