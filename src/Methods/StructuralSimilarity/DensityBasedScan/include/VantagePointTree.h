#ifndef VANTAGEPOINTTREE_H
#define VANTAGEPOINTTREE_H

#include <Eigen/Core>
#include <vector>
#include <queue>
#include <limits>
#include <random>
#include <spdlog/spdlog.h>

namespace Clustering {

    template<typename Scalar, Scalar ( *distance )(
            const Eigen::Matrix<Scalar,Eigen::Dynamic,1> &,
            const Eigen::Matrix<Scalar,Eigen::Dynamic,1> &)>
    class VantagePointTree {
        using VectorType = Eigen::Matrix<Scalar,Eigen::Dynamic,1>;

    public:
        VantagePointTree(const VantagePointTree&) = delete;
        VantagePointTree& operator=(const VantagePointTree&) = delete;

    private:

        struct HeapItem {
            HeapItem(size_t idx, Scalar dist)
                    : idx(idx), dist(dist) {
            }

            size_t idx;
            Scalar dist;

            bool operator<(const HeapItem &o) const {
                return dist < o.dist;
            }
        };

    public:
        static const Eigen::Index FIRTS_NODE_IDX = 1;

        explicit VantagePointTree(const std::vector<VectorType>& data, Scalar similarityDistance = 1e-7)
        :
        data_(data),
        similarityDistance(similarityDistance),
        rootIndex_(FIRTS_NODE_IDX),
        nextIndex_(FIRTS_NODE_IDX),
        randomDevice_(),
        mersenneTwister_(randomDevice_()) {}

        ~VantagePointTree() = default;

        void create() {
            rootIndex_ = FIRTS_NODE_IDX;
            nextIndex_ = FIRTS_NODE_IDX;
            //TODO DELdataset_ = dataset;

            //TODO DELconst auto &d = dataset_->data();

            nodelist_.resize(data_.size() + 1);
            itemsIndex_.resize(data_.size());

            for (size_t i = 0; i < data_.size(); ++i) {
                itemsIndex_[i] = i;
            }

            rootIndex_ = buildFromPoints(0, data_.size());
        }

        void searchByDistance(const VectorType &target, Scalar t, std::vector<std::pair<size_t, Scalar>> &nlist) const {
            nlist.clear();

            //TODO DEl const auto &d = dataset_->data();

            searchByDistance(rootIndex_, target, nlist, t);
        }

        void searchByK(const VectorType &target, size_t k, std::vector<std::pair<size_t, Scalar>> &nlist,
                       bool excludeExactQ = false) const {
            nlist.clear();

            Scalar t = std::numeric_limits<Scalar>::max();

            std::priority_queue<HeapItem> heap;

            searchByK(rootIndex_, target, nlist, k, heap, t, excludeExactQ);

            while (!heap.empty()) {
                const auto &top = heap.top();
                nlist.push_back(std::make_pair(top.idx, top.dist));
                heap.pop();
            }
        }

    private:
        struct Node {
            Eigen::Index index;
            Scalar threshold;
            Eigen::Index left;
            Eigen::Index right;

            Node() : index(0),
                     threshold(0.),
                     left(0),
                     right(0) {}
        };

        const std::vector<VectorType>& data_;
        Scalar similarityDistance;
        Eigen::Index rootIndex_;
        Eigen::Index nextIndex_;
        std::vector<Node> nodelist_;

        std::random_device randomDevice_;
        std::mt19937 mersenneTwister_;

        //std::shared_ptr<Dataset<Scalar>> dataset_;

        std::vector<size_t> itemsIndex_;

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

        Eigen::Index buildFromPoints(Eigen::Index lower, Eigen::Index upper) {
            if (upper == lower) {
                return 0;
            }

            Eigen::Index currentNodeIndex = nextIndex_++;

            Node &node = nodelist_[currentNodeIndex];
            node.index = lower;

            if (upper - lower > 1) {
                auto i = size_t(Scalar(mersenneTwister_()) / std::mt19937::max() * (upper - lower - 1)) + lower;

                std::swap(itemsIndex_[lower], itemsIndex_[i]);

                long median = (upper + lower) / 2;

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
                       /*const std::vector<VectorType> &d,*/
                       std::priority_queue<HeapItem> &heap,
                       Scalar &t,
                       bool excludeExactQ) const {
            if (nodeIndex == 0)
                return;

            const Node &node = nodelist_[nodeIndex];

            const Scalar dist = distance(data_[itemsIndex_[node.index]], target);

            if (dist < t) {
                if (!(excludeExactQ && dist < similarityDistance)) {

                    if (heap.size() == k) heap.pop();

                    heap.push(HeapItem(itemsIndex_[node.index], dist));

                    if (heap.size() == k) t = heap.top().dist;
                }
                /** do nothing on similar points **/
            }

            if (node.left == 0 && node.right == 0) {
                return;
            }

            if (dist < node.threshold) {
                if (dist - t <= node.threshold) {
                    searchByK(node.left, target, neighborList, k, heap, t, excludeExactQ);
                }

                if (dist + t >= node.threshold) {
                    searchByK(node.right, target, neighborList, k, heap, t, excludeExactQ);
                }

            } else {
                if (dist + t >= node.threshold) {
                    searchByK(node.right, target, neighborList, k, heap, t, excludeExactQ);
                }

                if (dist - t <= node.threshold) {
                    searchByK(node.left, target, neighborList, k, heap, t, excludeExactQ);
                }
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

            if (dist < t) {
                neighborList.push_back(std::make_pair(itemsIndex_[node.index], dist));
            }

            if (node.left == 0 && node.right == 0) {
                return;
            }

            if (dist < node.threshold) {
                if (dist - t <= node.threshold) {
                    searchByDistance(node.left, target, neighborList, t);
                }

                if (dist + t >= node.threshold) {
                    searchByDistance(node.right, target, neighborList, t);
                }

            } else {
                if (dist + t >= node.threshold) {
                    searchByDistance(node.right, target, neighborList, t);
                }

                if (dist - t <= node.threshold) {
                    searchByDistance(node.left, target, neighborList, t);
                }
            }
        }
    };
}

#endif // VANTAGEPOINTTREE_H
