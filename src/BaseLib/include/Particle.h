//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCPP_PARTICLE_H
#define AMOLQCPP_PARTICLE_H

#include <Eigen/Core>
#include <iostream>

enum class Type{
    none = 0
};

class Particle {
public:
    Particle(const Eigen::Vector3d& position, int type = int(Type::none));

    Eigen::Vector3d position() const;
    void position(const Eigen::Vector3d& position);

    virtual std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, const Particle& p);

    static double distance(const Particle &p1, const Particle &p2);

protected:
    virtual int charge() const;

    Eigen::Vector3d position_;
    int type_;
};

#endif //AMOLQCPP_PARTICLE_H
