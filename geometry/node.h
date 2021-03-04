#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers
#ifndef NODE_H
#define NODE_H

#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <cmath>
#include <Eigen/Core>
//#include "edge.h"

#include <concurrent_vector.h>
#include <concurrent_unordered_map.h>

namespace icy { class Node; }

class icy::Node
{
public:
    Node();
    void Reset();

    int id;          // sequential number in a given floe
    bool pinned = false;
    bool selected = false;
    double area;        // area that the node "represents", for applying various forces

//    tbb::concurrent_vector<icy::Element*> adjacent_elems;

    Eigen::Vector2d x_initial;  // initial configuration
    Eigen::Vector2d xn, vn;     // position and velocity at step n
    Eigen::Vector2d xt;         // at step n+1

};

#endif // NODE_H
#endif // Q_MOC_RUN
