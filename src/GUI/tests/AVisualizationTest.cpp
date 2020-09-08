// Copyright (C) 2018-2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
*/

#include <gmock/gmock.h>
#include "Visualization.h"

using namespace testing;
using namespace Eigen;
using namespace Visualization;

class AVisualizationTest : public Test {
public:
    ElectronsVectorCollection optimizationPath;

    void SetUp() override {
        Vector3d pos1(1,2,3);
        Vector3d pos2(-1,0,3.5);
        Vector3d pos3(-7.3,0.5,9);
        auto alpha = Spin::alpha;
        auto beta = Spin::beta;

        ElectronsVector electronsVector;

        for (int i = 0; i < 1000; i++){
            electronsVector = ElectronsVector();
            electronsVector.append(Electron(alpha,pos1));
            electronsVector.append(Electron(alpha,pos2));
            electronsVector.append(Electron(beta,pos3));

            optimizationPath.append(electronsVector);

            pos1(i%3) += 0.1;
            pos2(i%3) += 0.1;
            pos3(i%3) += 0.1;
        }
    }
};


TEST_F(AVisualizationTest, shortenPath) {
    auto shortenedPath = Visualization::shortenPath(optimizationPath,10);

    ASSERT_EQ(optimizationPath.typesVector().asEigenVector(),
              shortenedPath.typesVector().asEigenVector());
    ASSERT_EQ(optimizationPath[0].positionsVector().asEigenVector(),
              shortenedPath[0].positionsVector().asEigenVector());

    ASSERT_EQ(10,shortenedPath.numberOfEntities());

    ASSERT_EQ(optimizationPath[-1].positionsVector().asEigenVector(),
              shortenedPath[-1].positionsVector().asEigenVector());
}
