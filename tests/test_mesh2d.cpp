#include <gtest/gtest.h>
#include "../src/mesh2d.hpp"

TEST(Mesh2DTest, AddNode) {
    Mesh2D mesh;
    mesh.addNode(1.0, 2.0);
    mesh.addNode(3.0, 4.0);
    EXPECT_EQ(mesh.getNodes().size(), 2);
}

TEST(Mesh2DTest, AddTriangle) {
    Mesh2D mesh;
    mesh.addNode(0.0, 0.0);
    mesh.addNode(1.0, 0.0);
    mesh.addNode(0.0, 1.0);
    mesh.addTriangle(0, 1, 2);
    EXPECT_EQ(mesh.getTriangles().size(), 1);
}


TEST(Mesh2DTest, GenerateFacets) {
    Mesh2D mesh;
    mesh.addNode(0.0, 0.0);
    mesh.addNode(1.0, 0.0);
    mesh.addNode(0.0, 1.0);
    mesh.addTriangle(0, 1, 2);
    mesh.generateFacets();
    EXPECT_GT(mesh.getFacets().size(), 0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
