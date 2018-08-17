//
// Created by emmarose88: See https://github.com/emmarose88/HungarianAlgorithm
//

#ifndef LOCALSPINMULTIPLICITY_HUNGARIAN_H
#define LOCALSPINMULTIPLICITY_HUNGARIAN_H
#include <Eigen/Dense>

enum matchtype {MATCH_MIN=0,MATCH_MAX=1};

namespace Hungarian {
    void reduce(Eigen::MatrixXd& m);
    int hasMark(Eigen::VectorXd& v);
    void swapStarsAndPrimes(int i, int j, Eigen::MatrixXd& stars, Eigen::MatrixXd& primes);
    void findMatching(Eigen::MatrixXd &m, Eigen::MatrixXd & result, matchtype type);
}

#endif //LOCALSPINMULTIPLICITY_HUNGARIAN_H
