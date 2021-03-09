//#include <cmath>
//#include <cfloat>
//#include <algorithm>
#include "node.h"
//#include "model.h"
//#include "parameters_sim.h"
//#include "element.h"

icy::Node::Node()
{
    Reset();
}

void icy::Node::Reset()
{
    x_initial = xn = vn = xt = Eigen::Vector2d::Zero();
    id = -1;
    area = 0;

    //    adjacent_elems.reserve(7);
    //    adjacent_elems.clear();
}

void icy::Node::ComputeEquationEntries(EquationOfMotionSolver &eq, SimParams &prms, double timeStep)
{
    if(this->eqId<0) return;

    double mass = area * prms.Density * prms.Thickness;

    // for nodes that are not pinned, add the lumped mass matrix to the quadratic term of the equation
    Eigen::Matrix2d M_nd = Eigen::Matrix2d::Identity()*mass;
    eq.AddToQ(eqId, eqId, M_nd);

    Eigen::Vector2d lambda_n = xt-x_hat;
    Eigen::Vector2d linear_term = M_nd*lambda_n;
    eq.AddToC(eqId, linear_term);
}

