#include "mesh.h"

typedef Eigen::Vector2f vec2;
/*
    https://web.archive.org/web/20130126163405/http://geomalgorithms.com/a03-_inclusion.html
    > 1. an upward edge includes its starting endpoint, and excludes its final endpoint;
    > 2. a downward edge excludes its starting endpoint, and includes its final endpoint;
    > 3. horizontal edges are excluded
    > 4. the edge-ray intersection point must be strictly right of the point P.
    calculates winding order using those rules
    if a ray originating from a test point going in the +x direction hits a Edge
    slanting downwards it subtracts one from the winding number, and an upwards
    Edge adds one to the WN
*/
bool testPointInPoly(const vec2& P, const std::vector<Edge>& edges) {
    int wn = 0;
    for (const Edge& edge : edges) {
        const vec2& a = edge.a;
        const vec2& b = edge.b;

        if (a(1) <= P(1)) {
            if (b(1) > P(1))
                if (isLeft(a, b, P) > 0)
                    ++wn;
        }
        else {
            if (b(1) <= P(1))
                if (isLeft(a, b, P) < 0)
                    --wn;
        }
    }
    return wn;
}

std::vector<Edge> pointVecToEdge(const std::vector<vec2>& points) {
    std::vector<Edge> edges;
    edges.reserve(points.size());
    for (size_t i = 0; i < points.size(); ++i)
        edges.push_back(Edge(points[(i - 1) % points.size()], points[i]));
    return edges;
}

CartesianMesh::CartesianMesh(const std::vector<Edge> edges) {
    nx = 20; ny = 20;
    Bound polyBound;
    for (const auto& edge : edges) {
        polyBound = Union(polyBound, Bound(edge.a));
        polyBound = Union(polyBound, Bound(edge.b));
    }

    vec2 diag = polyBound.max - polyBound.min;
    dx = diag(0) / (float)(nx);
    dy = diag(1) / (float)(ny);

    typedef std::tuple<bool, size_t, size_t, size_t> cellInfo;
    std::vector<cellInfo> cellResults;
    cellResults.resize(nx * ny);

    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < ny; ++j) {
            const vec2 point = polyBound.min + vec2(
                ((float)i + 0.5)*dx,
                ((float)j + 0.5)*dy
            );
            if (testPointInPoly(point, edges)) {
                centroids.push_back(point);
                faceIDs.resize(faceIDs.size() + 1);
                cellResults[j + i * nx] = std::make_tuple(true, centroids.size() - 1, i, j);
            }
        }
    }

    // index 
    size_t n = 0;
    // maps neighboring cells
    for (size_t i = 0; i < nx; ++i) {
        for(size_t j = 0; j < ny; ++j){
            size_t index = j + i * nx;
            if (!std::get<bool>(cellResults[index]))
                continue;


            // check west side
            if (i != 0)
                if (std::get<bool>(cellResults[index - nx]))
                    faceIDs[index].W = std::make_pair(true, index - nx);
            
            // check east side
            if (i != nx - 1)
                if (std::get<bool>(cellResults[index + nx]))
                    faceIDs[index].E = std::make_pair(true, index + nx);

            // check south side
            if (j != 0)
                if (std::get<bool>(cellResults[index - 1]))
                    faceIDs[index].S = std::make_pair(true, index - 1);

            // check north side
            if (j != ny - 1)
                if (std::get<bool>(cellResults[index + 1]))
                    faceIDs[index].N = std::make_pair(true, index + 1);

        }
    }

}