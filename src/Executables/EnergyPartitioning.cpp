//
// Created by Michael Heuer on 28.08.18.
//

#include <RawDataReader.h>
#include "Sample.h"
#include "Logger.h"

int main(int argc, char *argv[]) {
    Logger::initialize();
    auto console = spdlog::get(Logger::name);

    std::set<Reference> references;
    std::vector<Sample> samples;

    RawDataReader reader(references,samples);
    reader.read("raw1.bin");

    for (auto ref : references) { ;
        console->info("{} {:03.6f}",ref.id_, ref.negLogSqrdProbabilityDensity_);
    }


}