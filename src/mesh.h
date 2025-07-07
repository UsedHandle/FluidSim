#pragma once

#include <vector>
#include <utility>
#include <Eigen/Dense>

#include "triangle.h"



constexpr float floatMax = std::numeric_limits<float>::max();
constexpr float floatMin = std::numeric_limits<float>::min();

struct Bound {
    typedef Eigen::Vector2f vec2;

    vec2 min, max;
    Bound() : min(floatMax, floatMax), max(floatMin, floatMin) {}
    Bound(vec2 a) : min(a), max(a) {}
    Bound(vec2 min, vec2 max) : min(min), max(max) {}
};

inline Bound Union(const Bound& a, const Bound& b) {
    typedef Eigen::Vector2f vec2;

    return Bound(vec2(std::min(a.min(0), b.min(0)),
                      std::min(a.min(1), b.min(1))),
                 vec2(std::max(a.max(0), b.max(0)),
                      std::max(a.max(1), b.max(1))));
}


static constexpr size_t North = 0;
static constexpr size_t West = 1;
static constexpr size_t East = 2;
static constexpr size_t South = 3;


struct faceID {
    static constexpr std::pair<bool, size_t> empty = std::make_pair(false, 0);
    union {
        struct {
            // first: true if N W E S side is another cell, false if it is dirchlet boundary condition
            // second: index of neighbor
            std::pair<bool, size_t> N, W, E, S;
        };
        std::array<std::pair<bool, size_t>, 4> faces;
    };
    faceID() : faces{ empty, empty, empty, empty } {}
};

struct CartesianMesh {
    std::vector<Eigen::Vector2f> centroids;
    std::vector<faceID> faceIDs;
    size_t nx, ny;
    float dx, dy;
    size_t known_start;

    CartesianMesh(const std::vector<Edge> edges);

    constexpr size_t size() const { return faceIDs.size(); }
};

// https://web.archive.org/web/20130126163405/http://geomalgorithms.com/a03-_inclusion.html
// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
inline float
isLeft(Eigen::Vector2f P0, Eigen::Vector2f P1, Eigen::Vector2f P2)
{
    return ((P1(0) - P0(0)) * (P2(1) - P0(1))
        - (P2(0) - P0(0)) * (P1(1) - P0(1)));
}

// tests if a point is in a polygon
// https://web.archive.org/web/20130126163405/http://geomalgorithms.com/a03-_inclusion.html
bool testPointInPoly(const Eigen::Vector2f& P, const std::vector<Edge>& edges);

// converts set of points to sets of edges  
std::vector<Edge> pointVecToEdge(const std::vector<Eigen::Vector2f>& points);