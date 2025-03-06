#include "mesh2d.hpp"
#include <iostream>
#include <set>

void Mesh2D::addNode(double x, double y) {
    nodes.push_back({x, y});
}

void Mesh2D::addTriangle(int n1, int n2, int n3) {
    triangles.push_back({{n1, n2, n3}});
}

void Mesh2D::generateFacets() {
    std::set<std::array<int, 2>> facetSet;

    for (int i = 0; i < triangles.size(); ++i) {
        std::array<std::array<int, 2>, 3> edges = {{
            {{triangles[i].nodes[0], triangles[i].nodes[1]}},
            {{triangles[i].nodes[1], triangles[i].nodes[2]}},
            {{triangles[i].nodes[2], triangles[i].nodes[0]}}
        }};

        for (auto& edge : edges) {
            if (edge[0] > edge[1]) std::swap(edge[0], edge[1]);  // Ensure ordering
            if (facetSet.insert(edge).second) {
                facets.push_back({static_cast<int>(facets.size()), edge, {i}, true});
            }
        }
    }
}

void Mesh2D::printMesh() const {
    std::cout << "Nodes:\n";
    for (size_t i = 0; i < nodes.size(); ++i)
        std::cout << i << ": (" << nodes[i].x << ", " << nodes[i].y << ")\n";

    std::cout << "\nTriangles:\n";
    for (size_t i = 0; i < triangles.size(); ++i)
        std::cout << i << ": " << triangles[i].nodes[0] << ", "
                  << triangles[i].nodes[1] << ", " << triangles[i].nodes[2] << "\n";
}
