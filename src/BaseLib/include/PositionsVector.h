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

    Interval(long start, long n)
    : start_(start),end_(start+n) {
        assert(start_ >= 0);
        assert(start_ <= end_);
    }
    explicit Interval(long idx)
    : Interval({idx,1}){}

    bool checkBounds(long numberOfEntities) const {
        return end() <= numberOfEntities;
    }

    long start() const { return start_;};
    long end() const { return end_;};
    long numberOfEntities() const { return end()-start();};

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

    PositionsRef positionsRef();
    void resetRef();

    void permute(long i, long j) override;
    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& permutation) override;
    void translate(const Eigen::Vector3d& shift);
    void rotateAroundOrigin(double angle, const Eigen::Vector3d &axisDirection);
    void rotate(double angle, const Eigen::Vector3d &center, const Eigen::Vector3d &axisDirection);

    PositionsVector(const PositionsVector& rhs);
    PositionsVector& operator=(const PositionsVector& rhs);

    PositionsVector& entity(long i);
    PositionsVector& slice(const Interval& interval);
    PositionsVector& all();

    friend std::ostream& operator<<(std::ostream& os, const PositionsVector& pc);

    const Eigen::VectorXd & positionsAsEigenVector() const;

    Eigen::VectorXd & positionsAsEigenVector();

private:
    Eigen::VectorXd positions_;
    const unsigned entityLength_ = 3;
    Interval sliceInterval_;
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
