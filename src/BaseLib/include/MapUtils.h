/* Copyright (C) 2020 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef INPSIGHTS_MAPUTILS_H
#define INPSIGHTS_MAPUTILS_H

#include <map>
#include <list>

namespace MapUtils {
    // finds all keys mapped to a value as a list of indices
    template<typename K, typename V>
    std::list<K> findByValue(std::map<K, V> mapOfElement, V value) {
        std::list<K> keyIndices;

        auto it = mapOfElement.begin();
        // Iterate through the map
        while (it != mapOfElement.end()) {
            // Check if value of this entry matches with given value
            if (it->second == value) {
                // Push the key in given map
                keyIndices.push_back(it->first);
            }
        // Go to next entry in map
        it++;
        }
        return keyIndices;
    }
}

#endif //INPSIGHTS_MAPUTILS_H
