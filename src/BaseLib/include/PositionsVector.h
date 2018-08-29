//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_POSITIONSVECTOR_H
#define AMOLQCPP_POSITIONSVECTOR_H

#include <Eigen/Core>
#include "AbstractVector.h"

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

    /* TODO
    void replace(long i);
    void remove(long i);
    ParticlesVector part(std::vector<long> indices);
    */

    friend std::ostream& operator<<(std::ostream& os, const PositionsVector& pc);

    const Eigen::VectorXd & positionsAsEigenVector() const;

    Eigen::VectorXd & positionsAsEigenVector();

    Eigen::Ref<Eigen::Vector3d> operator()(long i);

    const Eigen::Ref<const Eigen::Vector3d>& operator()(long i) const;

private:
    Eigen::VectorXd positions_;
    unsigned entityLength_ = 3;

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
