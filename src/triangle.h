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
bool pointInCircle(
    const vec2& P, 
    const vec2& A, 
    const vec2& B, 
    const vec2& C
);

struct QuartEdge {
    Eigen::Vector2f origin;
    QuartEdge* next;
    QuartEdge* rot;

    // returns the (dual) edge from the right to the left side
    // of the current edge
    inline QuartEdge* getRot() const { return rot; }

    // returns the edge with the same origin to the left of the current edge
    inline QuartEdge* getNext() const { return next; }

    // equivalent to two rotations
    // reverses direction of the edge
    inline QuartEdge* getSym() const { return rot->rot; }

    // returns the edge with the same origin to the right of the current edge
    inline QuartEdge* getPrev() const { return rot->next->rot;}

    // gets previous rotation
    inline QuartEdge* getPrevRot() const { return rot->rot->rot; }

    // Gets the next edge in counter clockwise order of polygon
    inline QuartEdge* getPolyNext() const { return getPrevRot()->next->rot; };

    inline const vec2& getOrigin() const { return origin; }
    inline vec2& getOrigin() { return origin; }

    inline const vec2& getDest() const { return getSym()->origin; }
    inline vec2& getDest() { return getSym()->origin; }

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

void swapNexts(QuartEdge* a, QuartEdge* b);

// places b's next as a's next and a's next as b's next
// additionally connects dual edges (the rot of a normal edge)
void splice(QuartEdge* a, QuartEdge* b);

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
void flip(QuartEdge* edge);
