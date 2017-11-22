//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_PARTICLE_H
#define AMOLQCGUI_PARTICLE_H

#include <Eigen/Core>

class Particle {
public:
    explicit Particle(const Eigen::Vector3d& position);
    Particle(double x, double y, double z);

    Eigen::Vector3d position() const;
    void position(const Eigen::Vector3d& position);

    static double distance(const Particle &p1, const Particle &p2);

protected:
    Eigen::Vector3d position_;
};
#endif //AMOLQCGUI_PARTICLE_H
