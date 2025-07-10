#pragma once

#include <array>

#include <Eigen/Dense>

typedef Eigen::Vector2f vec2;


struct Edge {
    vec2 a;
    vec2 b;

    Edge(vec2 a, vec2 b) : a(a), b(b) {}
};

/* Reference
*  https://ianthehenry.com/posts/delaunay/
*  Guibas & Stolfi "Primitives for the manipulation of general subdivisions and the computation of Voronoi" 1985
*/

// tests if point P is in the circumcircle formed by A, B and C
inline bool pointInCircle(
    const vec2& P,
    const vec2& A,
    const vec2& B,
    const vec2& C
)
{
    return Eigen::Matrix4f{
        {A(0), A(1), A(0) * A(0) + A(1) * A(1), 1},
        {B(0), B(1), B(0) * B(0) + B(1) * B(1), 1},
        {C(0), C(1), C(0) * C(0) + C(1) * C(1), 1},
        {P(0), P(1), P(0) * P(0) + P(1) * P(1), 1},
    }.determinant() > 0;
}

struct QuartEdge {
    Eigen::Vector2f origin;
    QuartEdge* next;
    QuartEdge* rot;

    // returns the (dual) edge from the right to the left side
    // of the current edge
    inline QuartEdge* getRot() { return rot; }
    inline const QuartEdge* getRot() const { return rot; }

    // returns the edge with the same origin to the left of the current edge
    inline QuartEdge* getNext() { return next; }
    inline const QuartEdge* getNext() const { return next; }

    // equivalent to two rotations
    // reverses direction of the edge
    inline QuartEdge* getSym() { return rot->rot; }
    inline const QuartEdge* getSym() const { return rot->rot; }

    // returns the edge with the same origin to the right of the current edge
    inline QuartEdge* getPrev() { return rot->next->rot;}
    inline const QuartEdge* getPrev() const { return rot->next->rot; }

    // gets previous rotation
    inline QuartEdge* getPrevRot() { return rot->rot->rot; }
    inline const QuartEdge* getPrevRot() const { return rot->rot->rot; }

    // Gets the edge that originates from this's dest and is connected to the face on this's left side
    inline QuartEdge* getPolyLNext() { return getPrevRot()->next->rot; };
    inline const QuartEdge* getPolyLNext() const { return getPrevRot()->next->rot; };

    // Gets the edge that originates from this's dest and is connected to the face on this's right side
    inline QuartEdge* getPolyRNext() { return rot->getPrev()->getPrevRot(); };
    inline const QuartEdge* getPolyRNext() const { return rot->getPrev()->getPrevRot(); };

    inline vec2& getOrigin() { return origin; }
    inline const vec2& getOrigin() const { return origin; }

    inline vec2& getDest() { return getSym()->origin; }
    inline const vec2& getDest() const { return getSym()->origin; }

    // initializes pointers with nullptr and origin with vec2(0.f, 0.f)
    QuartEdge() : next(nullptr), rot(nullptr), origin(0.f, 0.f) { }
    QuartEdge(
        vec2 start, vec2 end,
        QuartEdge* _rot = new QuartEdge(),
        QuartEdge* sym = new QuartEdge(),
        QuartEdge* prevRot = new QuartEdge()
    );
};

// contains original edge, rot, sym, and prevRot
struct QuadEdge {
    union {
        std::array<QuartEdge, 4> ref;
        struct {
            QuartEdge edge, rot, sym, prevRot;
        };
    };

    QuadEdge(vec2 start = vec2(0.f, 0.f), vec2 end = vec2(0.f, 0.f))
        : edge(QuartEdge(start, end, &rot, &sym, &prevRot)) {}
};

inline void
setStartEnd(QuartEdge* a, vec2 start, vec2 end) {
    a->origin = start;
    a->getDest() = end;
}

inline void swapNexts(QuartEdge* a, QuartEdge* b) {
    std::swap(a->next, b->next);
}

// places b's next as a's next and a's next as b's next
// additionally connects dual edges (the rot of a normal edge)
inline void splice(QuartEdge* a, QuartEdge* b) {
    swapNexts(a->next->rot, b->next->rot);
    swapNexts(a, b);
}

void makeTriangle(
    const vec2& A,
    const vec2& B,
    const vec2& C,
    QuadEdge* AB,
    QuadEdge* BC,
    QuadEdge* CA
);

// places edge between a and b and connects dual edges
void connect(QuartEdge* newEdge, QuartEdge* a, QuartEdge* b);


// removes edges connected to a and reverses connect
void sever(QuartEdge* a);

// flips diagonal edge of quadrilateral for delaunay triangulation
// the edge's origin is set to the point on its right (edge->getPrev()->getDest())
// the edge's dest is set to the point on its left (edge->getSym()->getPrev()->getDest())
void flip(QuartEdge* edge);

// places point in polygon and connects edges from the vertices of the polygon to the point
// the added edges are assigned to memory sequentially
// P must be to the left of polyEdge
void insertPoint(QuadEdge* memoryptr, QuartEdge* polyEdge, const vec2& P);


// returns positive number if left of, negative if right of, or 0 if on an segment from A to B
inline float pointLeftnessOfSegment(const vec2& A, const vec2& B, const vec2& P) {
    return Eigen::Matrix3f{
        {A(0), A(1), 1.f},
        {B(0), B(1), 1.f },
        {P(0), P(1), 1.f }
    }.determinant();
}

// returns if line segment AB intersects line segment CD once
// excludes colinearity
// excludes endpoints from intersection
inline bool lineSegmentIntersects(const vec2& A, const vec2& B, const vec2& C, const vec2& D) {
    // C and D are on opposite sides of the line formed by A and B,
    // A and B are on opposite sides of the line formed by C and D
    return (pointLeftnessOfSegment(A, B, C) > 0 != pointLeftnessOfSegment(A, B, D) > 0) &&
        (pointLeftnessOfSegment(C, D, A) > 0 != pointLeftnessOfSegment(C, D, B) > 0);
     
}
// returns positive number if left of, negative if right of, or 0 if on an edge
inline float pointLeftnessOfEdge(const QuartEdge* const edge, const vec2& P) {
    return pointLeftnessOfSegment(edge->origin, edge->getDest(), P);
}

// returns true if the edge is locally Delaunay
inline bool checkEdge(const QuartEdge* const edge) {
    const vec2& A = edge->getOrigin();
    const vec2& B = edge->getNext()->getDest();
    const vec2& C = edge->getDest();

    const vec2& P = edge->getPrev()->getDest();
    return !pointInCircle(P, A, B, C);
}

// checks if edge should be flipped
// bound1 is an "infinite" bounding triangle where its left side is the inside of the triangle
bool checkEdgeBoundaryAware(const QuartEdge* const edge, QuartEdge* bound1);

// will be stuck in infinite loop if the point is not bordering or inside a triangle
// returns false if point lies on an edge and true if search is successful
// starting edge is set to the edge of the triangle that contains the point,
// so that the edge's left side is the triangle's interior
template<typename QuartEdgePtr>
std::pair<bool, QuartEdgePtr> searchForTriangle(QuartEdgePtr startingEdge, const vec2& point) {
    bool isPrevEdgeChecked = false;
    auto* currentEdge = startingEdge;

    while (startingEdge != currentEdge || !isPrevEdgeChecked) {
        const float leftness = pointLeftnessOfEdge(currentEdge, point);
        if (leftness == 0.f) {
            return { false, currentEdge };
        }
        else if (leftness > 0.f) {
            currentEdge = currentEdge->getPolyLNext();
            isPrevEdgeChecked = true;
        }
        else if (leftness < 0.f) {
            currentEdge = currentEdge->getSym();
            startingEdge = currentEdge;
            isPrevEdgeChecked = false;
        }
    }
    return { true, currentEdge };
}