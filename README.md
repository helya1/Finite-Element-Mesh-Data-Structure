
![C++ CI](https://github.com/helya1/FiniteElementMesh/actions/workflows/ci.yml/badge.svg)

# Finite Element Mesh (2D)

A C++ implementation of a **2D finite element mesh** using triangular elements and explicit facets.

## Features
- Store **nodes** (points).
- Store **triangles** (elements).
- Generate **facets (edges)** with boundary detection.

## Build & Run
```bash
mkdir build
cd build
cmake ..
make
./FiniteElementMesh
