/* Copyright (C) 2017-2019 Michael Heuer.
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

#ifndef INPSIGHTS_POSITIONSVECTOR_H
#define INPSIGHTS_POSITIONSVECTOR_H

#include "InsertableVector.h"
#include <random>
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

    void shake(double radius, std::default_random_engine& rng);
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
