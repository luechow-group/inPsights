//
// Created by Morian Sonneton 16.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_PARTICLE_H
#define LOCALSPINMULTIPLICITY_PARTICLE_H
#include <Eigen/Dense>
#include <Eigen/Core>


class Particle {
public:
    Particle(double &x, double &y, double &z);
    virtual ~Particle();
    static double distance(const Particle &p1, const Particle &p2);
    const Eigen::Vector3d &getPosition() const;
    void setPosition(double &x, double &y, double &z);
private:
    Eigen::Vector3d position;
};


#endif //LOCALSPINMULTIPLICITY_PARTICLE_H
