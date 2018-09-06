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


int main(int argc, char *argv[]) {
    Logger::initialize();
    auto console = spdlog::get(Logger::name);

    std::set<Reference> references;
    std::vector<Sample> samples;

    RawDataReader reader(references,samples);
    reader.read("raw.bin");

    std::cout << "number of refs" << references.size() << std::endl;


    double increment, distThresh = 0.05;
    if(!references.empty()) {
        increment = (*references.rbegin()).negLogSqrdProbabilityDensity_ * 1e-4;
        console->info("increment: {}",increment);
        if(references.size() == 1){
            console->warn("No sorting because only one reference was found.");
            return true; // no sorting
        }
    } else {
        console->error("References are empty.");
        return false;
    }


    // initialize
    auto lit = references.begin();
    auto uit = references.begin();

    while (lit != references.end()){
        console->info("lit={}",(std::distance(references.begin(),lit)));
        uit = references.upper_bound(Reference((*lit).negLogSqrdProbabilityDensity_+increment));

        auto it = lit;
        it++;
        console->info("before inner while  it={}", (std::distance(references.begin(),it)));
        console->info("before inner while uit={}", (std::distance(references.begin(),uit)) );

        // check all refs
        while(it != uit) { // start with incremented iterator
            console->info("  it={}", (std::distance(references.begin(),it)) );

            // check hungarian
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
            console->info("end uit={}", (std::distance(references.begin(),uit)) );
        }
        console->info("after inner while  it={}", (std::distance(references.begin(),it)));
        console->info("after inner while uit={}", (std::distance(references.begin(),uit)) );
        lit = uit;
    }
    console->flush();

    for (auto a : references){
        std::cout << a.id_ << " {";
        for (auto b : a.associations_){
            std::cout << b << " ";
        }
        std::cout << " }" << std::endl;
    }

}