// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_MAPUTILS_H
#define INPSIGHTS_MAPUTILS_H

#include <map>
#include <list>

namespace MapUtils {
    // finds all keys mapped to a value as a set of indices
    template<typename K, typename V>
    std::set<K> findByValue(std::map<K, V> mapOfElement, V value) {
        std::set<K> keyIndices;

        auto it = mapOfElement.begin();
        // Iterate through the map
        while (it != mapOfElement.end()) {
            // Check if value of this entry matches with given value
            if (it->second == value) {
                // emplace the key in given map
                keyIndices.emplace(it->first);
            }
        // Go to next entry in map
        it++;
        }
        return keyIndices;
    }
}

#endif //INPSIGHTS_MAPUTILS_H
