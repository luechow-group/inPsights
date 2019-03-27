//
// Created by Michael Heuer on 2019-01-08.
//

#include <gmock/gmock.h>
#include <SurfaceData.h>
#include <iostream>
#include <yaml-cpp/yaml.h>
#include <EigenYamlConversion.h>

TEST(ASurfaceDataTest, YAMLConversion){
    SurfaceData data;
    data.vertices = {
            Vertex({1,2,3}),
            Vertex({4,5,6})};
    data.vertices[0].normal = {-1,-1,-1};

    data.triangles = {
            Triangle({0,1,2}),
            Triangle({2,3,0})
    };

    std::stringstream ss;
    ss << YAML::convert<SurfaceData>::encode(data);

    YAML::Emitter out;
    out << data;

    auto readData1 = YAML::Load(ss.str().c_str()).as<SurfaceData>();
    auto readData2 = YAML::Load(out.c_str()).as<SurfaceData>();

    ASSERT_EQ(data.vertices[0].position, readData1.vertices[0].position);
    ASSERT_EQ(data.vertices[1].position, readData1.vertices[1].position);
    ASSERT_EQ(data.vertices[0].position, readData2.vertices[0].position);
    ASSERT_EQ(data.vertices[1].position, readData2.vertices[1].position);

    ASSERT_EQ(data.vertices[0].normal, readData1.vertices[0].normal);
    ASSERT_EQ(data.vertices[1].normal, readData1.vertices[1].normal);
    ASSERT_EQ(data.vertices[0].normal, readData2.vertices[0].normal);
    ASSERT_EQ(data.vertices[1].normal, readData2.vertices[1].normal);

    ASSERT_EQ(data.triangles[0].indices, readData1.triangles[0].indices);
    ASSERT_EQ(data.triangles[1].indices, readData1.triangles[1].indices);
    ASSERT_EQ(data.triangles[0].indices, readData2.triangles[0].indices);
    ASSERT_EQ(data.triangles[1].indices, readData2.triangles[1].indices);

}