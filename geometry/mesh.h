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

    // nodes, elements and edges
    std::vector<icy::Node> nodes;
    std::vector<int> deformable_boundary_nodes;
    std::vector<icy::Element> elems;
    std::vector<std::pair<int,int>> boundary;

    std::vector<icy::Node> nodes_indenter;
    std::vector<std::pair<int,int>> boundary_indenter;

    // interaction with the indenter; populated in DetectContactPairs()
    struct Interaction
    {
        unsigned ndA_idx, ndB_idx, ndP_idx;
        double t, dist;
    };
    std::vector<Interaction> indenter_boundary_vs_deformable_nodes;
    std::vector<Interaction> deformable_boundary_vs_indenter_nodes;

    // at the "setting up the scene" stage - remesh 2d floe if needed
    void Reset(double CharacteristicLengthMax);
    void DetectContactPairs(double distance_threshold);
private:
    void GenrateIndenter(double CharacteristicLengthMax);
    double static SegmentPointDistance(Eigen::Vector2d A, Eigen::Vector2d B, Eigen::Vector2d P, double &t);
};
#endif
#endif // Q_MOC_RUN
