//
// Created by Michael Heuer on 29.10.17.
//

#ifndef AMOLQCGUI_PARTICLECOLLECTION_H
#define AMOLQCGUI_PARTICLECOLLECTION_H

#include "AbstractCollection.h"
#include "Particle.h"

class ParticleCollection : public AbstractCollection{
public:
    ParticleCollection();
    explicit ParticleCollection(const Eigen::VectorXd& positions);

    Particle operator[](long i) const;

    unsigned long numberOfParticles() const;

    void insert(const Particle& particle, long i);
    void append(const Particle& particle);
    void prepend(const Particle& particle);
    virtual void permute(long i, long j);

    /* TODO
    void replace(long i);
    void remove(long i);
    ParticleCollection part(std::vector<long> indices);
    */

    friend std::ostream& operator<<(std::ostream& os, const ParticleCollection& pc);

    Eigen::VectorXd positionsAsEigenVector() const;

protected:
    Eigen::VectorXd positions_;

private:
    long calculateStartIndex(long i) const;
};

#endif //AMOLQCGUI_PARTICLECOLLECTION_H
