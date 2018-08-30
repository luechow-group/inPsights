//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_POSITIONSVECTOR_H
#define AMOLQCPP_POSITIONSVECTOR_H

#include <Eigen/Core>
#include <memory>
#include "AbstractVector.h"
#include "Interval.h"
#include "ReturnAndReset.h"

using PositionsRef = Eigen::Ref<Eigen::VectorXd>;

class PositionsVector : public AbstractVector{
public:
    PositionsVector();
    explicit PositionsVector(const Eigen::VectorXd& positions);

    Eigen::Vector3d operator[](long i) const;

    PositionsVector& entity(long i, const Reset& resetType = Reset::Automatic);
    PositionsVector& slice(const Interval& interval, const Reset& resetType = Reset::Automatic);

    void insert(const Eigen::Vector3d& position, long i);
    void append(const Eigen::Vector3d& position);
    void prepend(const Eigen::Vector3d& position);
    void permute(long i, long j) override;

    PositionsRef positionsRef(const Usage& usage = Usage::NotFinished);
    Eigen::Vector3d position(long i, const Usage& usage = Usage::NotFinished);
    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& permutation, const Usage& usage);
    void permute(const Eigen::PermutationMatrix<Eigen::Dynamic>& permutation) override;
    void translate(const Eigen::Vector3d& shift, const Usage& usage = Usage::NotFinished);
    void rotateAroundOrigin(double angle, const Eigen::Vector3d &axisDirection, const Usage& usage = Usage::NotFinished);
    void rotate(double angle, const Eigen::Vector3d &center, const Eigen::Vector3d &axisDirection, const Usage& usage = Usage::NotFinished);

    PositionsVector(const PositionsVector& rhs);
    PositionsVector& operator=(const PositionsVector& rhs);

    const Eigen::VectorXd & positionsAsEigenVector() const;
    Eigen::VectorXd & positionsAsEigenVector();

    friend std::ostream& operator<<(std::ostream& os, const PositionsVector& pc);
    void resetRef();
    void resetter(const Usage& usage);

private:
    Eigen::VectorXd positions_;
    const unsigned entityLength_ = 3;
    Reset resetType_;
    Interval sliceInterval_;
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
