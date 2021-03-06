// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <Conversion.h>
#include <Eigen/Geometry>

std::vector<Vertex> Conversion::convertVertices(const std::vector<dualmc::Vertex>& dualMcVertices){
        std::vector<Vertex> vertices;

        for( const auto& v : dualMcVertices)
            vertices.emplace_back<Vertex>({{v.x, v.y, v.z}});

        return vertices;
    }

    std::pair<Triangle, Triangle> Conversion::quadToTrianglePair(const dualmc::Quad &quad) {
        return {Triangle({quad.i0, quad.i1, quad.i2}), Triangle({quad.i2, quad.i3, quad.i0})};
    }

    std::vector<Triangle> Conversion::quadsToTriangles(const std::vector<dualmc::Quad> &quads) {
        std::vector<Triangle> triangles;
        for (const auto &quad : quads) {
            auto trianglePair = quadToTrianglePair(quad);
            triangles.push_back(trianglePair.first);
            triangles.push_back(trianglePair.second);
        }

        return triangles;
    }

    Eigen::Vector3f Conversion::calculateFaceNormal(
            const Vertex &v1,
            const Vertex &v2,
            const Vertex &v3) {
        return (v2.position - v1.position).cross(v3.position - v1.position).normalized();
    }

    void Conversion::calculateVertexNormals(
            std::vector<Vertex> &vertices,
            const std::vector<Triangle> &triangles) {

        assert(vertices.size() > 3);
        assert(triangles.size() > 0);

        // add normals to vertices
        for (auto it = vertices.begin(); it != vertices.end(); ++it) {
            auto idx = static_cast<dualmc::QuadIndexType>(std::distance(vertices.begin(), it));
            Eigen::Vector3f vertexNormal = {0, 0, 0};

            // find triangles that contain the vertex
            for (const auto &t : triangles) {
                if (t.indices[0] == idx || t.indices[1] == idx || t.indices[2] == idx)
                    vertexNormal += calculateFaceNormal(
                            vertices[t.indices[0]],
                            vertices[t.indices[1]],
                            vertices[t.indices[2]]);
            }
            it.base()->normal = vertexNormal.normalized();
        }
    }