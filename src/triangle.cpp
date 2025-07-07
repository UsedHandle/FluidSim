#include "triangle.h"

#include <algorithm>

bool pointInCircle(
    const vec2& P,
    const vec2& A,
    const vec2& B,
    const vec2& C
)
{
    using namespace Eigen;
    const Matrix4f testMat{
        {A(0), A(1), A(0) * A(0) + A(1) * A(1), 1},
        {B(0), B(1), B(0) * B(0) + B(1) * B(1), 1},
        {C(0), C(1), C(0) * C(0) + C(1) * C(1), 1},
        {P(0), P(1), P(0) * P(0) + P(1) * P(1), 1},
    };
    return testMat.determinant() > 0;
}

QuartEdge::QuartEdge(
    vec2 start, vec2 end,
    QuartEdge* _rot,
    QuartEdge* sym,
    QuartEdge* prevRot
) : origin(start), rot(_rot) {
    sym->origin = end;
    rot->rot = sym;
    sym->rot = prevRot;
    prevRot->rot = this;

    // edge is not connected to anything
    next = this;
    sym->next = sym;

    // dual edges (cross the normal edge)
    // are part of the same face when there
    // is only one edge
    rot->next = prevRot;
    prevRot->next = rot;
}

void swapNexts(QuartEdge* a, QuartEdge* b) {
    std::swap(a->next, b->next);
}

void splice(QuartEdge* a, QuartEdge* b) {
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
 ){
    setStartEnd(&AB->edge, A, B);
    setStartEnd(&BC->edge, B, C);
    setStartEnd(&CA->edge, C, A);

    splice(AB->edge.getSym(), &BC->edge);
    splice(BC->edge.getSym(), &CA->edge);
    splice(CA->edge.getSym(), &AB->edge);
}

void connect(QuartEdge* newEdge, QuartEdge* a, QuartEdge* b) {
    setStartEnd(newEdge, a->getDest(), b->origin);
    splice(newEdge, a->getPolyNext());
    splice(newEdge->getSym(), b);
}

void sever(QuartEdge* a) {
    splice(a, a->getPrev());
    splice(a->getSym(), a->getSym()->getPrev());
}
#include<print>
void flip(QuartEdge* edge) {
    auto a = edge->getPrev();
    auto b = edge->getSym()->getPrev();
    std::print("{} {} {} {} {}", edge->getSym()->getOrigin(), a->getOrigin(), a->getDest(), b->getOrigin(), b->getDest());
    
    // severs edge from a and b->getSym()
    sever(edge);

    // flips edge
    connect(edge, a, b->getPolyNext());

    // updates origin and dest
    setStartEnd(edge, a->getDest(), b->getDest());
}
