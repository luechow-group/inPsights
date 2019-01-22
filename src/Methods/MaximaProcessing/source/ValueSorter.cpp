//
// Created by heuer on 19.12.18.
//

#include "ValueSorter.h"
#include <Logger.h>

bool ValueSorter::sortReferencesByValue(std::vector<Reference> &references) {
    using namespace Logger;
    if (references.empty()) {
        console->error("References are empty.");
        return false;
    } else if (references.size() == 1) {
        console->warn("No sorting because only one reference was found.");
        return true;
    }
    console->info("Sort references according to -ln|Ψ|² value...");
    std::sort(references.begin(), references.end());
    console->info("Finished value sorting.");
    return true;
}
