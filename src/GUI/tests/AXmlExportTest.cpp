//
// Created by heuer on 07.10.19.
//


#include <gmock/gmock.h>
#include "Cylinder.h"
#include "Sphere.h"
#include "X3dConverter.h"

using namespace testing;

class AXmlExportTest : public Test {
public:
    void SetUp() override {
    };
};

TEST_F(AXmlExportTest, out) {

    X3dConverter x3Dconverter("test.html");

    Qt3DCore::QEntity *root;

    QVector3D
        a = {0,-1,1},
        b = {1,2,1};

    Cylinder cylinder1(root, Qt::green,
                      {a,b},
                      0.2,
                      0.5);

    Sphere s1(root, Qt::red, a, 0.2, 0.7);
    Sphere s2(root, Qt::blue, b, 0.2, 0.3);

    x3Dconverter.addCylinder(cylinder1);
    x3Dconverter.addSphere(s1);
    x3Dconverter.addSphere(s2);

    x3Dconverter.closeScene();
}

