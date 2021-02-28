#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers

#ifndef FL333_H
#define FL333_H

#include <gmsh.h>

#include <vector>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <tuple>
#include <queue>
#include <set>

#include <QObject>
#include <QString>
#include <QDebug>

#include "element.h"

#include <concurrent_unordered_map.h>

namespace icy { class Mesh; class Node; class Element; class Edge;}

class icy::Mesh
{
public:
    Mesh(double CharacteristicLengthMax);

    // nodes, elements and edges
    std::vector<icy::Node> nodes;
    std::vector<icy::Element> elems;

    // at the "setting up the scene" stage - remesh 2d floe if needed
    void Reset(double CharacteristicLengthMax);

};
#endif
#endif // Q_MOC_RUN
