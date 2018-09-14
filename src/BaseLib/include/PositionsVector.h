//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_POSITIONSVECTOR_H
#define AMOLQCPP_POSITIONSVECTOR_H

#include "InsertableVector.h"
#include <Eigen/Core>


class PositionsVector : public InsertableVector<double>{
public:
    PositionsVector();
    explicit PositionsVector(const Eigen::VectorXd& positions);

    PositionsVector& entity(long i, const Reset& resetType = Reset::Automatic);
    PositionsVector& slice(const Interval& interval, const Reset& resetType = Reset::Automatic);

    void insert(const Eigen::Vector3d& position, long i);
    void append(const Eigen::Vector3d& position);
    void prepend(const Eigen::Vector3d& position);

    Eigen::Vector3d position(long i, const Usage& usage = Usage::Standard);
    void translate(const Eigen::Vector3d& shift, const Usage& usage = Usage::Standard);
    void rotateAroundOrigin(double angle, const Eigen::Vector3d &axisDirection, const Usage& usage = Usage::Standard);
    void rotate(double angle, const Eigen::Vector3d &center, const Eigen::Vector3d &axisDirection, const Usage& usage = Usage::Standard);


    Eigen::Vector3d operator[](long i) const;

    bool operator==(const PositionsVector& other) const;

    bool operator!=(const PositionsVector& other) const;

    friend std::ostream& operator<<(std::ostream& os, const PositionsVector& pc);
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
