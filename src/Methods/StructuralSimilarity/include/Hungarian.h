//
// Created by emmarose88: See https://github.com/emmarose88/HungarianAlgorithm
//
/*
BSD 3-Clause License

Copyright (c) 2017, Emma Rose Streshinsky
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
//
// Modified by Michael Heuer on 08.04.18 based on http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
//


#ifndef AMOLQCPP_HUNGARIAN_H
#define AMOLQCPP_HUNGARIAN_H

#include <Eigen/Dense>
#include <limits>

enum class Matchtype {MIN, MAX};

template<typename Scalar>
class Hungarian {
public:
    static Eigen::PermutationMatrix<Eigen::Dynamic>
    findMatching(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &m, Matchtype type = Matchtype::MIN) {
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> n = m; // make a copy of m for reducing
        Eigen::Index dim = n.rows(); // dimension of matrix, used for checking if we've reduced the matrix enough yet

        Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> stars(m.rows(), m.cols()); // matrix for storing our "starred" 0s (0*)
        stars.fill(false);
        Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> primes(m.rows(), m.cols()); // matrix for storing our "primed" 0s (0')
        primes.fill(false);
        Eigen::Array<bool, Eigen::Dynamic, 1> rowCover(m.rows()); // keep track of which rows are "covered"
        rowCover.fill(false);
        Eigen::Array<bool, Eigen::Dynamic, 1> colCover(m.cols()); // keep track of which columns are "covered"
        colCover.fill(false);

        // to do maximization rather than minimization, we have to
        // transform the matrix by subtracting every value from the maximum
        if (type == Matchtype::MAX) {
            Scalar max = n.maxCoeff();
            Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> maxMat(n.rows(), n.cols());
            maxMat.fill(max);
            n = maxMat - n;
        }

        // Step 1
        // Reduce matrix
        reduce(n);

        // Step 2
        // Find a zero in the matrix. If there is no starred zero in
        // its row or column, star Z. Repeat for each element in the matrix.
        for (int i = 0; i < n.rows(); i++) {
            for (int j = 0; j < n.cols(); j++) {
                if (n(i, j) == 0 && !rowCover(i) && !colCover(j)) {
                    stars(i, j) = true;
                    rowCover(i) = true;
                    colCover(j) = true;
                }
            }
        }
        // covers need to be cleared for following steps
        rowCover.fill(false);
        colCover.fill(false);

        while (true) {
            // Step 3
            // Cover all columns that have a starred zero
            // If the number of columns with starred zeroes equals the matrix
            // dimensions, we are done! Otherwise, move on to step 4.
            step3:
            for (int j = 0; j < n.cols(); j++) {
                Eigen::Array<bool, Eigen::Dynamic, 1> col = stars.col(j);
                if (hasMark(col) >= 0) {
                    colCover(j) = true;
                }
            }

            Eigen::Index coveredCols = 0;
            for (int i = 0; i < colCover.size(); ++i) {
                if (colCover[i]) coveredCols++;
            }

            if (coveredCols == dim) {
                Eigen::VectorXi permutation(dim);
                for (Eigen::Index i = 0; i < dim; i++) {
                    for (Eigen::Index j = 0; j < dim; j++) {
                        if (stars(i, j)) {
                            permutation[i] = j;
                        }
                    }
                }
                return Eigen::PermutationMatrix<Eigen::Dynamic>(permutation);
            }

            // Step 4
            // Find a non-covered zero and prime it
            step4:
            for (int i = 0; i < n.rows(); i++) {
                for (int j = 0; j < n.cols(); j++) {
                    if (n(i, j) == 0 && !rowCover(i) && !colCover(j)) {
                        primes(i, j) = true;
                        // if no starred zero in the row...
                        Eigen::Array<bool, Eigen::Dynamic, 1> row = stars.row(i);
                        if (hasMark(row) < 0) {
                            // Step 5
                            // swap stars and primes
                            swapStarsAndPrimes(i, j, stars, primes);

                            // clear lines
                            rowCover.fill(false);
                            colCover.fill(false);

                            goto step3;
                        } else {
                            // cover row
                            rowCover(i) = true;

                            // uncover column of the starred zero in the same row
                            int col = hasMark(row);
                            colCover(col) = false;
                        }
                    }
                }
            }

            // Step 6
            // Should now be no more uncovered zeroes
            // Get the minimum uncovered element
            Scalar min = std::numeric_limits<Scalar>::max();
            for (int i = 0; i < n.rows(); i++) {
                for (int j = 0; j < n.cols(); j++) {
                    if (!rowCover(i) && !colCover(j) && n(i, j) < min) {
                        min = n(i, j);
                    }
                }
            }

            // Subtract minimum from uncovered elements, add it to elements covered twice
            for (int i = 0; i < n.rows(); i++) {
                for (int j = 0; j < n.cols(); j++) {
                    if (!rowCover(i) && !colCover(j)) {
                        n(i, j) -= min;
                    } else if (rowCover(i) && colCover(j)) {
                        n(i, j) += min;
                    }
                }
            }

            goto step4;
        }
    };

private:
    static void reduce(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> &m) {
        // subtract row minimum from each row
        for (int i = 0; i < m.rows(); i++) {
            Scalar minElement = m.row(i).minCoeff();
            Eigen::Matrix<Scalar, Eigen::Dynamic, 1> rMinusMin(m.rows());
            rMinusMin.fill(-minElement);
            m.row(i) += rMinusMin;
        }
    }

    static int hasMark(Eigen::Array<bool, Eigen::Dynamic, 1> &v) {
        for (int i = 0; i < v.size(); i++) {
            if (v(i)) {
                return i;
            }
        }
        return -1;
    };

    static void swapStarsAndPrimes(int i, int j, Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> &stars, Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> &primes) {
        int primeRow = i;
        int primeCol = j;

        bool done = false;
        while (!done) {
            // find row index of row that has a 0* in the same col as the current 0'
            Eigen::Array<bool, Eigen::Dynamic, 1> col = stars.col(primeCol);
            int starInPrimeColRow = hasMark(col);

            if (starInPrimeColRow < 0) {
                // star the prime we're looking at
                primes(primeRow, primeCol) = false;
                stars(primeRow, primeCol) = true;
                done = true;
            } else {
                // find which col has a 0' in the same row as z1
                Eigen::Array<bool, Eigen::Dynamic, 1> row = primes.row(starInPrimeColRow);
                int primeInStarRowCol = hasMark(row);

                // star first primed zero
                primes(primeRow, primeCol) = false;
                stars(primeRow, primeCol) = true;

                // unstar starred zero
                stars(starInPrimeColRow, primeCol) = false;

                // set index of last prime, will check it's column for 0*s next
                primeRow = starInPrimeColRow;
                primeCol = primeInStarRowCol;
            }
        }
        // clear primes
        primes.fill(false);
    };
};

#endif //AMOLQCPP_HUNGARIAN_H
