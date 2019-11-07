/* Copyright (C) 2019 Michael Heuer.
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

