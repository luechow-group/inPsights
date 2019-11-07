/* Copyright (C) 2019 Michael Heuer.
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

#ifndef INPSIGHTS_LINKEDPARTICLE_H
#define INPSIGHTS_LINKEDPARTICLE_H

#include <Eigen/Core>

template <typename Type> class LinkedParticle{
public:
    LinkedParticle(
            const Eigen::Ref<Eigen::Vector3d> &positionRef,
            int* typeNameRef)
            : position_(positionRef), type_(typeNameRef)
    {}

    Eigen::Ref<Eigen::Vector3d>& positionRef(){
        return position_;
    }

    Eigen::Vector3d position() const {
        return Eigen::Vector3d(position_);
    }

    virtual void setPosition(const Eigen::Vector3d & position){
        position_ = position;
    }

    Type type() const {
        return static_cast<Type>(*type_);
    }

    int* typeRef(){
        return type_;
    }

    void setType(const Type & type){
        *type_ = static_cast<int>(type);
    }

private:
    Eigen::Ref<Eigen::Vector3d> position_;
    int* type_;
};

#endif //INPSIGHTS_LINKEDPARTICLE_H
