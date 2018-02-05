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

    friend std::ostream& operator<< (std::ostream& os, const Particle& p);

    static double distance(const Particle &p1, const Particle &p2);

protected:
    Eigen::Vector3d position_;
};

namespace ParticleFormat{
    static const unsigned significantDigits = 6;
    static const std::string separator = " ";
    static const Eigen::IOFormat particleFormat = Eigen::IOFormat(significantDigits, 0, separator, "\n", "", "", "", "");
};

#endif //AMOLQCGUI_PARTICLE_H
