#include "triangle.h"

#include <algorithm>



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
    splice(newEdge, a->getPolyLNext());
    splice(newEdge->getSym(), b);
}

void sever(QuartEdge* a) {
    splice(a, a->getPrev());
    splice(a->getSym(), a->getSym()->getPrev());
}

void flip(QuartEdge* edge) {
    auto a = edge->getPrev();
    auto b = edge->getSym()->getPrev();
    
    // severs edge from a and b->getSym()
    sever(edge);

    // flips edge
    connect(edge, a, b->getPolyLNext());

    // updates origin and dest
    setStartEnd(edge, a->getDest(), b->getDest());
}

void insertPoint(QuadEdge* memoryptr, QuartEdge* polyEdge, const vec2& P) {
    const auto firstSpoke = memoryptr;
    setStartEnd(&memoryptr->edge, polyEdge->getOrigin(), P);

    splice(&memoryptr->edge, polyEdge);
    do {
        const auto newSpoke = memoryptr+1;
        // memoryptr is the previous spoke created
        connect(&(newSpoke)->edge, polyEdge, memoryptr->edge.getSym());
        polyEdge = newSpoke->edge.getPrev();
        memoryptr = newSpoke;
    } while (polyEdge->getPolyLNext() != &firstSpoke->edge);
}

bool checkEdgeBoundaryAware(const QuartEdge* const edge, QuartEdge* bound1) {
    for (int i = 0; i < 3; ++i) {
        // checks if edge is a boundary edge
        if (edge == bound1 || edge == bound1->getSym())
            return false;

        // edge connects to boundary vertex
        if (edge->origin == bound1->origin || edge->getDest() == bound1->origin) {
            const bool flipIntersectsEdge = lineSegmentIntersects(
                edge->origin, edge->getDest(),
                edge->getNext()->getDest(), edge->getPrev()->getDest()
            );
            // if the flipped edge does not intersect the original, an inside out triangle will be created
            // when the edge is flipped
            if (!flipIntersectsEdge)
                return false;
            else return true;
        }

        // flipped edge connects to boundary vertex
        const bool flipTouchesBoundary = edge->getNext()->getDest() == bound1->origin ||
            edge->getSym()->getNext()->getDest() == bound1->origin;
        if (flipTouchesBoundary)
            return false;
        bound1 = bound1->getPolyRNext();
    }
    return checkEdge(edge);
}

