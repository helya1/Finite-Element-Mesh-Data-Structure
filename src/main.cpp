#include "mesh2d.hpp"

int main() {
    Mesh2D mesh;

    mesh.addNode(0.0, 0.0);
    mesh.addNode(1.0, 0.0);
    mesh.addNode(0.0, 1.0);

    mesh.addTriangle(0, 1, 2);

    mesh.generateFacets();
    mesh.printMesh();

    return 0;
}
