#include "element.h"
//#include "model.h"
#include "node.h"

#include <cstdlib>
#include <algorithm>

#include <iostream>
#include <vector>
#include <utility>

#include <QtDebug>


void icy::Element::Reset(void)
{
    for(int i=0;i<3;i++) {
        nds[i] = nullptr;
        adj_elems[i] = nullptr;
    }
    area_initial = 0;
}

void icy::Element::PrecomputeInitialArea()
{
    Eigen::Matrix2d J;
    J << nds[0]->x_initial.x()-nds[2]->x_initial.x(), nds[1]->x_initial.x()-nds[2]->x_initial.x(),
            nds[0]->x_initial.y()-nds[2]->x_initial.y(), nds[1]->x_initial.y()-nds[2]->x_initial.y();
    area_initial = J.determinant()/2;
}


void icy::Element::AddToSparsityStructure(EquationOfMotionSolver &eq)
{
    // register the positions of non-zero entries with the EquationOfMotionSolver
    eq.AddElementToStructure(nds[0]->eqId, nds[1]->eqId);
    eq.AddElementToStructure(nds[0]->eqId, nds[2]->eqId);
    eq.AddElementToStructure(nds[1]->eqId, nds[2]->eqId);
}


void icy::Element::ComputeEquationEntries(EquationOfMotionSolver &eq, SimParams &prms, double timeStep)
{
    // evaluate quadratic, linear and constant terms of the equation of motion

    // assemble the result into the EquationOfMotionSolver

    // for testing
    SpringModel(eq, prms, timeStep, nds[0],nds[1]);
    SpringModel(eq, prms, timeStep, nds[1],nds[2]);
    SpringModel(eq, prms, timeStep, nds[2],nds[0]);
}


void icy::Element::SpringModel(EquationOfMotionSolver &eq, SimParams &prms, double timeStep, Node *nd1, Node *nd2)
{
    double d = (nd1->x_initial - nd2->x_initial).norm();
    double d_new = (nd1->xt - nd2->xt).norm();
    double d_new_sq = d_new*d_new;
    double d_new_cube = d_new_sq * d_new;

    double x1 = nd1->xt.x();
    double x2 = nd2->xt.x();
    double y1 = nd1->xt.y();
    double y2 = nd2->xt.y();


    Eigen::Vector2d c1;
    c1(0) = x2-x1;
    c1(1) = y2-y1;
    c1*=timeStep*timeStep*prms.YoungsModulus*(d-d_new)/d_new;

    Eigen::Vector2d c2 = -c1;
    eq.AddToC(nd1->eqId, c1);
    eq.AddToC(nd2->eqId, c2);

    Eigen::Matrix2d Q11;
    Q11(0,0) = 1+d*(-d_new_sq+(x1-x2)*(x1-x2))/d_new_cube;
    Q11(0,1)=Q11(1,0)= d*(x1-x2)*(y1-y2)/d_new_cube;
    Q11(1,1) = 1+d*(-d_new_sq+(y1-y2)*(y1-y2))/d_new_cube;
    Q11*=timeStep*timeStep*prms.YoungsModulus;

    // Q22 === Q11
    // Q12 === -Q11

    Eigen::Matrix2d Q12 = -Q11;
    eq.AddToQ(nd1->eqId, nd1->eqId, Q11);
    eq.AddToQ(nd2->eqId, nd2->eqId, Q11);
    eq.AddToQ(nd1->eqId, nd2->eqId, Q12);
    eq.AddToQ(nd2->eqId, nd1->eqId, Q12);
}


/*


    // assemble
    for(int i=0;i<3;i++) {
        int row = nds[i]->lsId;
        Eigen::Matrix<double,DOFS,1> locF = F.block(i*DOFS,0,DOFS,1);
        ls.SubtractRHS(row, locF);
        for(int j=0;j<3;j++) {
            int col = nds[j]->lsId;
            Eigen::Matrix<double,DOFS,DOFS> loc_dF = dF.block(i*DOFS,j*DOFS,DOFS,DOFS);
            ls.AddLHS(row, col, loc_dF);
        }
    }
}
*/

