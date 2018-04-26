//
// Created by Leonard Reuter on 15.03.18.
//

#include <gtest/gtest.h>
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
        auto alpha = Spins::SpinType::alpha;
        auto beta = Spins::SpinType::beta;

        ElectronsVector electronsVector;

        for (int i = 0; i < 1000; i++){
            electronsVector = ElectronsVector();
            electronsVector.append(Electron(pos1,alpha));
            electronsVector.append(Electron(pos2,alpha));
            electronsVector.append(Electron(pos3,beta));

            optimizationPath.append(electronsVector);

            pos1(i%3) += 0.1;
            pos2(i%3) += 0.1;
            pos3(i%3) += 0.1;
        }
    }
};

TEST_F(AVisualizationTest, empty) {
}

TEST_F(AVisualizationTest, shortenPath) {
    auto shortenedPath = Visualization::shortenPath(optimizationPath,10);

    ASSERT_EQ(optimizationPath.spinTypesVector().typesAsEigenVector()(),
              shortenedPath.spinTypesVector().typesAsEigenVector()());
    ASSERT_EQ(optimizationPath[0].positionsVector().positionsAsEigenVector(),
              shortenedPath[0].positionsVector().positionsAsEigenVector());

    ASSERT_EQ(10,shortenedPath.numberOfEntities());

    ASSERT_EQ(optimizationPath[-1].positionsVector().positionsAsEigenVector(),
              shortenedPath[-1].positionsVector().positionsAsEigenVector());
}
