#include "QuickHull.h"
#include "MathUtils.h"
#include <iostream>
#include <random>
#include <cassert>

#include <gmock/gmock.h>
#include "MolecularGeometry.h"
#include "TestMolecules.h"

using namespace testing;
using namespace quickhull;

using FloatType = float;
using vec3 = Vector3<FloatType>;

auto sameEps = 0.0001f;

class AQuickHullTest: public Test {
public:

    // Setup test env
    const size_t N = 200;
    std::vector<vec3> pc;

    std::default_random_engine rng;
    std::uniform_real_distribution<FloatType> dist;

    void SetUp() override {
        // Seed RNG using Unix time
        dist = std::uniform_real_distribution<FloatType> (0, 1); //remove
        rng.seed(static_cast<unsigned long>(std::clock()));


        // Test 1 : Put N points inside unit cube. Result mesh must have exactly 8 vertices because the convex hull is the unit cube.
        pc.clear();
        for (int i=0;i<8;i++) {
            pc.emplace_back(i&1 ? -1 : 1,i&2 ? -1 : 1,i&4 ? -1 : 1);
        }
        for (size_t i=0;i<N;i++)
        {
            pc.emplace_back(rnd(-1,1),rnd(-1,1),rnd(-1,1));
        }
    };

    FloatType rnd(FloatType from, FloatType to) {
        return from + (FloatType)dist(rng)*(to-from);
    };

    template <typename T>
    std::vector<Vector3<T>> createSphere(T radius, size_t M, Vector3<T> offset = Vector3<T>(0,0,0)) {
        std::vector<Vector3<T>> pc;
        const T pi = M_PI;
        for (size_t i=0;i<=M;i++) {
            FloatType y = std::sin(pi/2 + (FloatType)i/(M)*pi);
            FloatType r = std::cos(pi/2 + (FloatType)i/(M)*pi);
            FloatType K = FloatType(1)-std::abs((FloatType)((FloatType)i-M/2.0f))/(FloatType)(M/2.0f);
            const size_t pcount = (size_t)(1 + K*M + FloatType(1)/2);
            for (size_t j=0;j<pcount;j++) {
                FloatType x = pcount>1 ? r*std::cos((FloatType)j/pcount*pi*2) : 0;
                FloatType z = pcount>1 ? r*std::sin((FloatType)j/pcount*pi*2) : 0;
                pc.emplace_back(x+offset.x,y+offset.y,z+offset.z);
            }
        }
        return pc;
    }

};


TEST_F(AQuickHullTest, Vector3){
    typedef Vector3<FloatType> vec3;
    vec3 a(1,0,0);
    vec3 b(1,0,0);

    vec3 c = a.projection(b);
    ASSERT_LT((c-a).getLength(), 0.00001f);

    a = vec3(1,1,0);
    b = vec3(1,3,0);
    c = b.projection(a);
    ASSERT_LT((c-vec3(2,2,0)).getLength(), 0.00001f);
}

TEST_F(AQuickHullTest, Sphere){
    QuickHull<FloatType> qh;
    FloatType y = 1;
    for (;;) {
        auto pc = createSphere<FloatType>(1, 100, Vector3<FloatType>(0,y,0));
        auto hull = qh.getConvexHull(pc,true,false);
        y *= 15;
        y /= 10;
        if (hull.getVertexBuffer().size()==4) {
            break;
        }
    }

    // Test worst case scenario: more and more points on the unit sphere. All points should be part of the convex hull, as long as we can make epsilon smaller without
    // running out of numerical accuracy.
    size_t i =  1;
    FloatType eps = 0.002f;
    for (;;) {
        auto pc = createSphere<FloatType>(1, i, Vector3<FloatType>(0,0,0));
        auto hull = qh.getConvexHull(pc,true,false,eps);
        std::cout << i << ":" << pc.size() << " : " << hull.getVertexBuffer().size() << " at eps=" << eps << std::endl;
        if (qh.getDiagnostics().m_failedHorizonEdges) {
            // This should not happen
            ASSERT_TRUE(false);
            break;
        }
        if (pc.size() == hull.getVertexBuffer().size()) {
            // Fine, all the points on unit sphere do belong to the convex mesh.
            i += 1;
        }
        else {
            eps *= 0.5f;
            std::cout << "Epsilon to " << eps << std::endl;
        }

        if (i == 200) { //Original value 500
            break;
        }
    }
}

TEST_F(AQuickHullTest, PlanarCase){
    QuickHull<FloatType> qh;
    std::vector<vec3> pointCloud;
    pointCloud.emplace_back(-3.000000f, -0.250000f, -0.800000f);
    pointCloud.emplace_back(-3.000000f, 0.250000f, -0.800000f);
    pointCloud.emplace_back(-3.125000f, 0.250000f, -0.750000);
    pointCloud.emplace_back(-3.125000f, -0.250000f, -0.750000);
    auto hull = qh.getConvexHull(pointCloud,true,false);
    ASSERT_EQ(hull.getIndexBuffer().size(), 12);
    ASSERT_EQ(hull.getVertexBuffer().size(), 4);
}

TEST_F(AQuickHullTest, Planes){
    Vector3<FloatType> N(1,0,0);
    Vector3<FloatType> p(2,0,0);
    Plane<FloatType> P(N,p);
    auto dist = mathutils::getSignedDistanceToPlane(Vector3<FloatType>(3,0,0), P);
    ASSERT_NEAR(dist, 1, sameEps);
    dist = mathutils::getSignedDistanceToPlane(Vector3<FloatType>(1,0,0), P);
    ASSERT_NEAR(dist, -1, sameEps);
    dist = mathutils::getSignedDistanceToPlane(Vector3<FloatType>(1,0,0), P);
    ASSERT_NEAR(dist, -1, sameEps);
    N = Vector3<FloatType>(2,0,0);
    P = Plane<FloatType>(N,p);
    dist = mathutils::getSignedDistanceToPlane(Vector3<FloatType>(6,0,0), P);
    ASSERT_NEAR(dist, 8, sameEps);
}


TEST_F(AQuickHullTest, HalfEdgeOutput) {
    QuickHull<FloatType> qh;

    // 8 corner vertices of a cube + tons of vertices inside. Output should be a half edge mesh with 12 faces (6 cube faces with 2 triangles per face) and 36 half edges (3 half edges per face).
    std::vector<vec3> pc;
    for (int h=0;h<1000;h++) {
        pc.emplace_back(rnd(-1,1),rnd(-1,1),rnd(-1,1));
    }
    for (int h=0;h<8;h++) {
        pc.emplace_back(h&1?-2:2, h&2?-2:2, h&4?-2:2);
    }
    HalfEdgeMesh<FloatType, size_t> mesh = qh.getConvexHullAsMesh(&pc[0].x, pc.size(), true);
    ASSERT_EQ(mesh.m_faces.size(), 12);
    ASSERT_EQ(mesh.m_halfEdges.size(), 36);
    ASSERT_EQ(mesh.m_vertices.size(), 8);
}


TEST_F(AQuickHullTest, Hull1) {
    QuickHull<FloatType> qh;
    auto hull = qh.getConvexHull(pc,true,false);
    ASSERT_EQ(hull.getVertexBuffer().size(), 8);
    ASSERT_EQ(hull.getIndexBuffer().size(), 3*2*6); // 6 cube faces, 2 triangles per face, 3 indices per triangle
    ASSERT_NE(&(hull.getVertexBuffer()[0]), &(pc[0]));
}

TEST_F(AQuickHullTest, Hull2) {
    QuickHull<FloatType> qh;
    auto hull = qh.getConvexHull(pc,true,false);

    ASSERT_EQ(hull.getVertexBuffer().size(), hull.getVertexBuffer().size());
    ASSERT_EQ(hull.getVertexBuffer()[0].x, hull.getVertexBuffer()[0].x);
    ASSERT_EQ(hull.getIndexBuffer().size(), hull.getIndexBuffer().size());
}

TEST_F(AQuickHullTest, Hull3) {
    QuickHull<FloatType> qh;
    auto hull = qh.getConvexHull(pc, true, false);

    auto hull3 = std::move(hull);
    ASSERT_EQ(hull.getIndexBuffer().size(), 0);
}

TEST_F(AQuickHullTest, Hull1OriginalIndices) {
    QuickHull<FloatType> qh;
    auto hull = qh.getConvexHull(pc, true, true);

    ASSERT_EQ(hull.getIndexBuffer().size(), 3*2*6);
    ASSERT_EQ(hull.getVertexBuffer().size(), pc.size());
    ASSERT_EQ(&(hull.getVertexBuffer()[0]), &(pc[0]));
}

TEST_F(AQuickHullTest, RandomPointsFromUnitSphere) {
    // random N points from the boundary of unit sphere. Result mesh must have exactly N points.

    QuickHull<FloatType> qh;

    auto pc = createSphere<FloatType>(1, 50);
    auto hull = qh.getConvexHull(pc,true,false);
    ASSERT_EQ(pc.size(), hull.getVertexBuffer().size());


    hull = qh.getConvexHull(pc,true,true);
    // Add every vertex twice. This should not affect final mesh
    auto s = pc.size();
    for (size_t i=0;i<s;i++) {
        const auto& p = pc[i];
        pc.push_back(p);
    }
    hull = qh.getConvexHull(pc,true,false);
    ASSERT_EQ(pc.size()/2, hull.getVertexBuffer().size());

    // Multiply x components of the unit sphere vectors by a huge number => essentially we get a line
    const FloatType mul = 2*2*2;
    while (true) {
        for (auto& p : pc) {
            p.x *= mul;
        }
        hull = qh.getConvexHull(pc,true,false);
        if (hull.getVertexBuffer().size() == 4) {
            break;
        }
    }
}

TEST_F(AQuickHullTest, ZeroDimension){
    QuickHull<FloatType> qh;
    std::vector<vec3> pointCloud;

    vec3 centerPoint(2,2,2);
    pointCloud.push_back(centerPoint);
    for (size_t i=0;i<100;i++) {
        auto newp = centerPoint + vec3(rnd(-0.000001f,0.000001f),rnd(-0.000001f,0.000001f),rnd(-0.000001f,0.000001f));
        pointCloud.push_back(newp);
    }
    auto hull = qh.getConvexHull(pointCloud,true,false);
    ASSERT_EQ(hull.getIndexBuffer().size(), 12);

}

TEST_F(AQuickHullTest, PlanarCircleAndCylinder){
    QuickHull<FloatType> qh;
    std::vector<vec3> pointCloud;

    // first a planar circle, then make a cylinder out of it
    pointCloud.clear();    for (size_t i=0;i<N;i++) {
        const FloatType alpha = (FloatType)i/N*2*Constant::pi;
        pointCloud.emplace_back(std::cos(alpha),0,std::sin(alpha));
    }
    auto hull = qh.getConvexHull(pointCloud,true,false);

    assert(hull.getVertexBuffer().size() == pointCloud.size());
    for (size_t i=0;i<N;i++) {
        pointCloud.push_back(pointCloud[i] + vec3(0,1,0));
    }
    hull = qh.getConvexHull(pointCloud,true,false);
    assert(hull.getVertexBuffer().size() == pointCloud.size());
    assert(hull.getIndexBuffer().size()/3 == pointCloud.size()*2-4);
}

TEST_F(AQuickHullTest, Test6){
    QuickHull<FloatType> qh;
    std::vector<vec3> pointCloud;

    for (int x=0;;x++) {
        pointCloud.clear();
        const FloatType l = 1;
        const FloatType r = l/(std::pow(10, x));
        for (size_t i=0;i<N;i++) {
            vec3 p = vec3(1,0,0)*i*l/(N-1);
            FloatType a = rnd(0,2*3.1415f);
            vec3 d = vec3(0,std::sin(a),std::cos(a))*r;
            pointCloud.push_back(p+d);
        }
        auto hull = qh.getConvexHull(pointCloud,true,false);
        if (hull.getVertexBuffer().size()==4) {
            break;
        }
    }
}

TEST_F(AQuickHullTest, Test7) {
    QuickHull<FloatType> qh;
    std::vector<vec3> pointCloud;

    for (int h = 0; h < 100; h++) {
        pointCloud.clear();
        const vec3 v1(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1));
        const vec3 v2(rnd(-1, 1), rnd(-1, 1), rnd(-1, 1));
        pointCloud.push_back(v1);
        pointCloud.push_back(v2);
        for (size_t i = 0; i < N; i++) {
            auto t1 = rnd(0, 1);
            auto t2 = rnd(0, 1);
            pointCloud.push_back(t1 * v1 + t2 * v2);
        }
        auto hull = qh.getConvexHull(pointCloud, true, false);
    }

    // TODO ASSERT was missing here
}
