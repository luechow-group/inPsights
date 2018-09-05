//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include "Sample.h"
#include "Logger.h"


void reduce(){

    // in equal range
    // do hungarian
    // if ok, unite sets and average centroid
}

int main(int argc, char *argv[]) {
    Logger::initialize();
    auto console = spdlog::get(Logger::name);

    std::set<Reference> references;
    std::vector<Sample> samples;

    RawDataReader reader(references,samples);
    reader.read("raw.bin");

    //for (auto ref : references) { ;
    //    console->info("{} {:03.6f}",ref.id_, ref.negLogSqrdProbabilityDensity_);
    //}



    double increment;
    if(!references.empty()) {
        increment = (*references.rbegin()).negLogSqrdProbabilityDensity_ * 1e-5;
    } else {
        console->error("References are empty.");
        return false;
    }
    console->info("increment: {}",increment);
    console->flush();


    // first range
    auto uit = references.begin();
    while (uit != references.end()){

        auto lit = uit;
        uit = references.upper_bound(Reference((*lit).negLogSqrdProbabilityDensity_+increment));


        console->info("{} ------ {}",(*lit).id_, (*lit).negLogSqrdProbabilityDensity_);

        for(auto it = lit; it != uit; it++){
            console->info("{}", (*it).id_);
        }

        console->info("{} ------ {}",(*uit).id_, (*uit).negLogSqrdProbabilityDensity_);

    }
    console->flush();

    //auto lit = references.lower_bound(Reference(435.1)); // fake ref for comparison
    //auto uit = references.upper_bound(Reference(435.3)); // fake ref for comparison


}