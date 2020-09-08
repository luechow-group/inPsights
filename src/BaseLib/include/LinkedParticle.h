// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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
