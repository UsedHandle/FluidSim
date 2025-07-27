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

/*
QuartEdge::QuartEdge(QuartEdge&& other) noexcept :
    origin(std::move(other.origin)),
    next(std::move(other.next)),
    rot(std::move(other.rot))
{
    if(getPrev() != &other)
        getPrev()->next = this;
    else {
        rot->next->rot = this;
        next = this;
    }
    getPrevRot()->rot = this;
}

QuartEdge& QuartEdge::operator=(QuartEdge&& other) noexcept {
    if (this != &other) {
        origin = std::move(other.origin);
        next = std::move(other.next);
        rot = std::move(other.rot);
        getPrev()->next = this;
        getPrevRot()->rot = this;

        other.origin = vec2(nanf("0"), nanf("0"));
        other.rot = nullptr;
        other.next = nullptr;
    }

    return *this;
}*/
QuadEdge::QuadEdge(QuadEdge&& other) noexcept
    : edge(QuartEdge(other.edge.origin, other.sym.origin, &rot, &sym, &prevRot)){
    splice(&edge, &other.edge);
    splice(other.edge.getPrev(), &other.edge);
    splice(&sym, &other.sym);
    splice(other.sym.getPrev(), &other.sym);
}

QuadEdge& QuadEdge::operator=(QuadEdge&& other) noexcept {
    if (this != &other){
        edge.origin = std::move(other.edge.origin);
        sym.origin = std::move(other.sym.origin);
        splice(&edge, &other.edge);
        splice(other.edge.getPrev(), &other.edge);
        splice(&sym, &other.sym);
        splice(other.sym.getPrev(), &other.sym);
    }
    
    return *this;
}

QuadEdge::~QuadEdge() {
    destroy(&edge);
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

size_t insertPoint(QuadEdge* memoryptr, QuartEdge* polyEdge, const vec2& P) {
    const auto firstSpoke = memoryptr;
    setStartEnd(&memoryptr->edge, polyEdge->getOrigin(), P);

    size_t edgesCreated = 1;
    splice(&memoryptr->edge, polyEdge);
    do {
        const auto newSpoke = memoryptr+1;
        // memoryptr is the previous spoke created
        connect(&(newSpoke)->edge, polyEdge, memoryptr->edge.getSym());
        polyEdge = newSpoke->edge.getPrev();
        memoryptr = newSpoke;
        ++edgesCreated;
    } while (polyEdge->getPolyLNext() != &firstSpoke->edge);
    return edgesCreated;
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

std::pair<size_t, QuartEdge*> insertPointDelaunay(QuadEdge* memoryPtr, QuartEdge* searchStart, QuartEdge* bound1, const vec2& point) {
    QuartEdge* polyBound;

    auto results = searchForTriangle(searchStart, point);
    if (!std::get<bool>(results)) {
        polyBound = results.second->getPrev();
        destroy(results.second);
    }
    else {
        polyBound = results.second;
    }
    size_t edgesCreated = insertPoint(memoryPtr, polyBound, point);

    QuartEdge* firstEdge = polyBound;
    // makes sure that firstEdge is unflippable
    while (checkEdgeBoundaryAware(firstEdge, bound1)) {
        flip(firstEdge);
        firstEdge = firstEdge->getNext()->getSym();
    }
    polyBound = firstEdge->getPolyLNext()->getPrev();
    // loops around ring
    while (polyBound != firstEdge) {
        if (checkEdgeBoundaryAware(polyBound, bound1)) {
            flip(polyBound);
            polyBound = polyBound->getNext()->getSym();
            continue;
        }
        polyBound = polyBound->getPolyLNext()->getPrev();
    }

    return { edgesCreated, firstEdge };
}

bool isPointInBoundaryTriangle(const vec2& point, QuartEdge* bound1) {
    for (int i = 0; i < 3; ++i) {
        if (pointLeftnessOfEdge(bound1, point) > 0)
            continue;
        return false;
    }
    return true;
}

std::array<vec2, 3> getTriPoints(QuartEdge* bound) {
    return std::array<vec2, 3>{
        bound->getOrigin(),
        bound->getPolyLNext()->getOrigin(),
        bound->getPolyLNext()->getPolyLNext()->getOrigin()
    };
}