#if !defined(Q_MOC_RUN) // MOC has a glitch when parsing TBB headers
#ifndef NODE_H
#define NODE_H

#include <concurrent_vector.h>
//#include <concurrent_unordered_map.h>

//#include <vector>
//#include <algorithm>
//#include <map>
//#include <unordered_map>
//#include <cmath>
#include <Eigen/Core>

#include "equationofmotionsolver.h"
#include "parameters_sim.h"

namespace icy { class Node; }

class icy::Node
{
public:
    Node();
    void Reset();

    int id, eqId;       // sequential number of a node; identificator in the equation of motion (if not pinned)
    bool pinned = false;
    bool selected = false;
    double area;        // area that the node "represents", for applying various forces

//    tbb::concurrent_vector<icy::Element*> adjacent_elems;

    Eigen::Vector2d x_initial;  // initial configuration
    Eigen::Vector2d xn, vn;     // position and velocity at step n
    Eigen::Vector2d xt;         // tentative coordinates
    Eigen::Vector2d x_hat;

    void ComputeEquationEntries(EquationOfMotionSolver &eq, SimParams &prms, double timeStep);

};

#endif // NODE_H
#endif // Q_MOC_RUN
