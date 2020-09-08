// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later


#include <gmock/gmock.h>
#include "Cylinder.h"
#include "Sphere.h"
#include "X3domConverter.h"

using namespace testing;

class AXmlExportTest : public Test {
public:
    void SetUp() override {
    };
};

TEST_F(AXmlExportTest, XmlTestFile) {

    X3domConverter x3Dconverter("test.html");
    QVector3D
        a = {0,-1,1},
        b = {1,2,1};

    Cylinder cylinder1(nullptr, Qt::green,
                      {a,b},
                      0.2,
                      0.5);

    Sphere s1(nullptr, Qt::red, a, 0.2, 0.7);
    Sphere s2(nullptr, Qt::blue, b, 0.2, 0.3);

    x3Dconverter.addCylinder(cylinder1);
    x3Dconverter.addSphere(s1);
    x3Dconverter.addSphere(s2);

    x3Dconverter.closeScene();
}

