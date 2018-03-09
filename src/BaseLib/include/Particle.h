//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_PARTICLE_H
#define AMOLQCGUI_PARTICLE_H

#include <Eigen/Core>
#include <iostream>

class Particle {
public:
    explicit Particle(const Eigen::Vector3d& position);
    Particle(double x, double y, double z);

    Eigen::Vector3d position() const;
    void position(const Eigen::Vector3d& position);

    virtual std::string toString() const;
    friend std::ostream& operator<< (std::ostream& os, const Particle& p);

    static double distance(const Particle &p1, const Particle &p2);

protected:
    virtual int charge() const;

    Eigen::Vector3d position_;
};

#endif //AMOLQCGUI_PARTICLE_H
