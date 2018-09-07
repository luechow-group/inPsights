//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include <Sample.h>
#include <Logger.h>
#include <HungarianHelper.h>


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

        /*PUT INTO METHOD*/
        //TODO CHECK MULTIPLICITY
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
    double distThresh = 0.01;

    std::set<Reference> references;
    std::vector<Sample> samples;

    RawDataReader reader(references,samples);
    reader.read("raw.bin");

    console->info("number of refs {}",references.size());

    GlobalIdentiySorter globalIdentiySorter(references, samples, distThresh);
    globalIdentiySorter.doSort();

    console->flush();

    /*for (auto a : references){
        std::cout << a.id_ << " {";
        for (auto b : a.associations_){
            std::cout << b << " ";
        }
        std::cout << " }" << std::endl;
    }*/


    std::set<SimilarReferencesCollection> similarReferencesCollections;
    // add first ref

    for(auto& r : references){

        // check for similarity

        bool isSimilarQ = false;
        for(auto& sr : similarReferencesCollections){
            // check if it belongs to some


            // calc perm r->sr so that we can store them in sr
            auto costMatrix = Metrics::positionalDistances(
                    r.maximum_.positionsVector(),
                    sr.representativeReference_->maximum_.positionsVector());
            auto bestMatch = Hungarian<double>::findMatching(costMatrix);
            auto dist = mostDeviatingParticleDistance(
                    r.maximum_.positionsVector(),
                    bestMatch,
                    sr.representativeReference_->maximum_.positionsVector());

            if(dist < distThresh) {
                auto a = SimilarReference(std::make_unique<Reference>(r),bestMatch);
                sr.similarReferences_.emplace(a);
                isSimilarQ = true;
            }
        }



        if(!isSimilarQ) similarReferencesCollections.emplace(r);
    }



}

/*
 *
std::vector<Type> v = ....;
std::string myString = ....;
auto it = find_if(v.begin(), v.end(), [&myString](const Type& obj) {return obj.getName() == myString;})

if (it != v.end())
{
  // found element. it is an iterator to the first matching element.
  // if you really need the index, you can also get it:
  auto index = std::distance(v.begin(), it);
}

 */
