#ifndef MESH2D_HPP
#define MESH2D_HPP

#include <vector>
#include <array>

struct Facet {
    int id;                        // Unique identifier
    std::array<int, 2> nodes;       // Nodes defining the facet
    std::vector<int> elements;      // Elements sharing this facet
    bool isBoundary;               // Boundary flag
};

struct Triangle {
    std::array<int, 3> nodes;       // Triangle has 3 nodes
};

struct Node {
    double x, y;                   // Node coordinates
};

class Mesh2D {
private:
    std::vector<Node> nodes;
    std::vector<Triangle> triangles;
    std::vector<Facet> facets;

public:
    Mesh2D() = default;

    void addNode(double x, double y);
    void addTriangle(int n1, int n2, int n3);
    void generateFacets();
    void printMesh() const;
};

#endif
