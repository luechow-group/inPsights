//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_POSITIONSVECTOR_H
#define AMOLQCPP_POSITIONSVECTOR_H

#include <Eigen/Core>
#include <memory>
#include "AbstractVector.h"

using PositionsRef = Eigen::Ref<Eigen::VectorXd>;

class Interval {
public:
    Interval() = default;

    Interval(long start, long end)
    : start_(start),end_(end) {
        assert(start_ >= 0);
        assert(start_ <= end_);
    }
    explicit Interval(long element)
    : Interval({element,element}){}

    bool checkBounds(long numberOfEntities) const {
        return end() <= numberOfEntities;
    }

    long start() const { return start_;};
    long end() const { return end_;};
    long numberOfEntities() const { return end()-start()+1;}; //TODO empty intervals?

private:
    long start_,end_;
};

class PositionsVector : public AbstractVector{
public:
    PositionsVector();
    explicit PositionsVector(const Eigen::VectorXd& positions);

    Eigen::Vector3d operator[](long i) const;

    void insert(const Eigen::Vector3d& position, long i);
    void append(const Eigen::Vector3d& position);
    void prepend(const Eigen::Vector3d& position);
    void permute(long i, long j) override;
    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& permutation) override;

    PositionsVector(const PositionsVector& rhs);
    PositionsVector& operator=(const PositionsVector& rhs);

    PositionsVector& slice(long i);
    PositionsVector& slice(const Interval& interval);
    PositionsVector& all();

    friend std::ostream& operator<<(std::ostream& os, const PositionsVector& pc);

    const Eigen::VectorXd & positionsAsEigenVector() const;

    Eigen::VectorXd & positionsAsEigenVector();

    Eigen::Ref<Eigen::Vector3d> operator()(long i);//TODO DEPRECATED

    const Eigen::Ref<const Eigen::Vector3d>& operator()(long i) const;

private:
    Eigen::VectorXd positions_;
    const unsigned entityLength_ = 3;
public:
    std::unique_ptr<PositionsRef> positionsRefPtr_;

    long calculateIndex(long i) const override ;

    Eigen::PermutationMatrix<Eigen::Dynamic> adaptedToEntityLength(const Eigen::PermutationMatrix<Eigen::Dynamic> &permutation);

};

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<> struct convert<PositionsVector> {
        static Node encode(const PositionsVector &rhs);
        static bool decode(const Node &node, PositionsVector &rhs);
    };
    Emitter &operator<<(Emitter &out, const PositionsVector &p) ;
}


#endif //AMOLQCPP_POSITIONSVECTOR_H
