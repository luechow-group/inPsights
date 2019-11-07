/*
 * Copyright (c) 2017, Vladimir Skvortsov
 * All rights reserved.
 *
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of the FreeBSD Project.
 */
/*
 * Copyright (C) 2018 Michael Heuer.
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

#ifndef VANTAGEPOINTTREE_H
#define VANTAGEPOINTTREE_H

#include <Eigen/Core>
#include <vector>
#include <queue>
#include <limits>
#include <random>
#include <spdlog/spdlog.h>


template<typename Scalar,
        typename VectorType,
        Scalar (*distance)(const VectorType &,const VectorType &)>
class VantagePointTree {

    struct HeapItem {
        HeapItem(size_t idx, Scalar dist)
                : idx(idx), dist(dist) {}

        size_t idx;
        Scalar dist;

        bool operator<(const HeapItem &o) const {
            return dist < o.dist;
        }
    };

    struct DistanceComparator {
        size_t itemIndex;
        const std::vector<VectorType> &items;

        DistanceComparator(size_t i, const std::vector<VectorType> &data)
                : itemIndex(i), items(data) {
        }

        bool operator()(size_t a, size_t b) {
            return distance(items[itemIndex], items[a]) < distance(items[itemIndex], items[b]);
        }
    };

    struct Node {
        Eigen::Index index, left, right;
        Scalar threshold;

        Node()
        : index(0), left(0), right(0), threshold(0.) {}
    };

public:
    VantagePointTree(const VantagePointTree &) = delete;

    VantagePointTree &operator=(const VantagePointTree &) = delete;

    explicit VantagePointTree(const std::vector<VectorType> &data, Scalar similarityDistance)
            :
            data_(data),
            similarityDistance(similarityDistance),
            rootIndex_(1),
            nextIndex_(1),
            randomDevice_(),
            mersenneTwister_(randomDevice_()) {
        nodelist_.resize(data_.size() + 1);
        itemsIndex_.resize(data_.size());

        for (size_t i = 0; i < data_.size(); ++i)
            itemsIndex_[i] = i;

        rootIndex_ = buildFromPoints(0, data_.size());
    }

    ~VantagePointTree() = default;

    void searchByDistance(
            const VectorType &target, Scalar tau,
            std::vector<std::pair<size_t, Scalar>> &neighborList) const {

        neighborList.clear();
        searchByDistance(rootIndex_, target, neighborList, tau);
    }

    void searchByK(const VectorType &target, size_t k, std::vector<std::pair<size_t, Scalar>> &neighborList,
                   bool excludeExactQ = false) const {
        neighborList.clear();

        Scalar tau = std::numeric_limits<Scalar>::max();

        std::priority_queue<HeapItem> heap;

        searchByK(rootIndex_, target, neighborList, k, heap, tau, excludeExactQ);

        while (!heap.empty()) {
            const auto &top = heap.top();
            neighborList.push_back(std::make_pair(top.idx, top.dist));
            heap.pop();
        }
    }

private:
    const std::vector<VectorType> &data_;
    Scalar similarityDistance;
    Eigen::Index rootIndex_;
    Eigen::Index nextIndex_;
    std::vector<Node> nodelist_;

    std::random_device randomDevice_;
    std::mt19937 mersenneTwister_;

    std::vector<size_t> itemsIndex_;


    Eigen::Index buildFromPoints(Eigen::Index lower, Eigen::Index upper) {
        if (upper == lower)
            return 0;

        Eigen::Index currentNodeIndex = nextIndex_++;

        Node &node = nodelist_[currentNodeIndex];
        node.index = lower;

        if (upper - lower > 1) {
            auto i = size_t(Scalar(mersenneTwister_()) / std::mt19937::max() * (upper - lower - 1)) + lower;

            std::swap(itemsIndex_[lower], itemsIndex_[i]);

            Eigen::Index median = (upper + lower) / 2;

            std::nth_element(
                    itemsIndex_.begin() + lower + 1,
                    itemsIndex_.begin() + median,
                    itemsIndex_.begin() + upper,
                    DistanceComparator(itemsIndex_[lower], data_));
            node.threshold = distance(data_[itemsIndex_[lower]], data_[itemsIndex_[median]]);

            node.index = lower;
            node.left = buildFromPoints(lower + 1, median);
            node.right = buildFromPoints(median, upper);
        }

        return currentNodeIndex;
    }

    void searchByK(Eigen::Index nodeIndex,
                   const VectorType &target,
                   std::vector<std::pair<size_t, Scalar>> &neighborList,
                   size_t k,
                   std::priority_queue<HeapItem> &heap,
                   Scalar &tau,
                   bool excludeExactQ) const {

        if (nodeIndex == 0)
            return;

        const Node &node = nodelist_[nodeIndex];

        const Scalar dist = distance(data_[itemsIndex_[node.index]], target);

        if (dist < tau) {
            if (!(excludeExactQ && dist < similarityDistance)) {

                if (heap.size() == k) heap.pop();
                heap.push(HeapItem(itemsIndex_[node.index], dist));
                if (heap.size() == k) tau = heap.top().dist;
            }
            /* do nothing on similar points */
        }

        if (node.left == 0 && node.right == 0)
            return;

        if (dist < node.threshold) {
            if (dist - tau <= node.threshold)
                searchByK(node.left, target, neighborList, k, heap, tau, excludeExactQ);

            if (dist + tau >= node.threshold)
                searchByK(node.right, target, neighborList, k, heap, tau, excludeExactQ);

        } else {
            if (dist + tau >= node.threshold)
                searchByK(node.right, target, neighborList, k, heap, tau, excludeExactQ);

            if (dist - tau <= node.threshold)
                searchByK(node.left, target, neighborList, k, heap, tau, excludeExactQ);
        }
    }

    void searchByDistance(Eigen::Index nodeIndex,
                          const VectorType &target,
                          std::vector<std::pair<size_t, Scalar>> &neighborList,
                          Scalar t) const {
        if (nodeIndex == 0)
            return;

        const Node &node = nodelist_[nodeIndex];

        const Scalar dist = distance(data_[itemsIndex_[node.index]], target);

        if (dist < t)
            neighborList.push_back(std::make_pair(itemsIndex_[node.index], dist));

        if (node.left == 0 && node.right == 0)
            return;

        if (dist < node.threshold) {
            if (dist - t <= node.threshold)
                searchByDistance(node.left, target, neighborList, t);

            if (dist + t >= node.threshold)
                searchByDistance(node.right, target, neighborList, t);
        } else {
            if (dist + t >= node.threshold)
                searchByDistance(node.right, target, neighborList, t);

            if (dist - t <= node.threshold)
                searchByDistance(node.left, target, neighborList, t);
        }
    }
};


#endif // VANTAGEPOINTTREE_H
