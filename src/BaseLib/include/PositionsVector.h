//
// Created by Michael Heuer on 29.10.17.
//

#ifndef INPSIGHTS_POSITIONSVECTOR_H
#define INPSIGHTS_POSITIONSVECTOR_H

#include "InsertableVector.h"
#include <Eigen/Core>

class PositionsVector : public InsertableVector<double>{
public:
    PositionsVector();
    explicit PositionsVector(const Eigen::VectorXd& positions);

    void insert(const Eigen::Vector3d& position, long i);
    void append(const Eigen::Vector3d& position);
    void prepend(const Eigen::Vector3d& position);

    Eigen::Vector3d position(long i);

    void translate(const Eigen::Vector3d& shift);
    void rotateAroundOrigin(double angle, const Eigen::Vector3d &axisDirection);
    void rotate(double angle, const Eigen::Vector3d &center, const Eigen::Vector3d &axisDirection);

    Eigen::Vector3d operator[](long i) const;

    bool operator==(const PositionsVector& other) const;

    bool operator!=(const PositionsVector& other) const;

    friend std::ostream& operator<<(std::ostream& os, const PositionsVector& pc);

    void shake(double radius);
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


#endif //INPSIGHTS_POSITIONSVECTOR_H
