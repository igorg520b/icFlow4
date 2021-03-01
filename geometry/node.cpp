#include <cmath>
#include <cfloat>
#include <algorithm>
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


