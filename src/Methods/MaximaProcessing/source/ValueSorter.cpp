//
// Created by heuer on 19.12.18.
//

#include "ValueSorter.h"
#include <spdlog/spdlog.h>

bool ValueSorter::sortReferencesByValue(std::vector<Reference> &references) {
    if (references.empty()) {
        spdlog::error("References are empty.");
        return false;
    } else if (references.size() == 1) {
        spdlog::warn("No sorting because only one reference was found.");
        return true;
    }
    spdlog::info("Sort references according to -ln|Ψ|² value...");
    std::sort(references.begin(), references.end());
    spdlog::info("Finished value sorting.");
    return true;
}
