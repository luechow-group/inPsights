//
// Created by Moria on 16.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_PARTICLE_H
#define LOCALSPINMULTIPLICITY_PARTICLE_H
#include <vector>


class Particle {
public:
    Particle(double x, double y, double z);
    virtual ~Particle();
    static double distance(const Particle &p1, const Particle &p2);
    const std::vector<double> &get_position() const;
private:
    std::vector<double> position;
};


#endif //LOCALSPINMULTIPLICITY_PARTICLE_H
