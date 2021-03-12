#include "element.h"
#include "node.h"

#include <cstdlib>
#include <algorithm>

#include <iostream>
#include <vector>
#include <utility>
#include <cmath>

// #include <QtDebug>


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
    NeoHookeanElasticity(eq, prms, timeStep);

    // for testing
//    SpringModel(eq, prms, timeStep, nds[0],nds[1]);
//    SpringModel(eq, prms, timeStep, nds[1],nds[2]);
//    SpringModel(eq, prms, timeStep, nds[2],nds[0]);
}

void icy::Element::NeoHookeanElasticity(EquationOfMotionSolver &eq, SimParams &prms, double h)
{
    double E = prms.YoungsModulus;
    double nu = prms.PoissonsRatio;
    double lambda = (E*nu)/((1.0+nu)*(1.0-2.0*nu));
    double K = E/(3.0*(1.0-2.0*nu));
    double mu = (K-lambda)*3.0/2.0;

    double X1 = nds[0]->x_initial.x();
    double X2 = nds[1]->x_initial.x();
    double X3 = nds[2]->x_initial.x();
    double Y1 = nds[0]->x_initial.y();
    double Y2 = nds[1]->x_initial.y();
    double Y3 = nds[2]->x_initial.y();

    double x1 = nds[0]->xt.x();
    double x2 = nds[1]->xt.x();
    double x3 = nds[2]->xt.x();
    double y1 = nds[0]->xt.y();
    double y2 = nds[1]->xt.y();
    double y3 = nds[2]->xt.y();

    Eigen::Matrix2d Dm, Dm_inv, Ds, F, Finv, FT, FinvT;
    Dm << X1-X3, X2-X3, Y1-Y3, Y2-Y3; // reference shape matrix
    double W = prms.Thickness*Dm.determinant()/2;
    Dm_inv = Dm.inverse();

    Ds << x1-x3, x2-x3, y1-y3, y2-y3; // deformed shape matrix
    F = Ds*Dm.inverse();    // deformation gradient
    FT = F.transpose();
    Finv = F.inverse();
    FinvT = Finv.transpose();

    // derivatives of Ds with respect to x1,y1,x2,y2,x3,y3
    Eigen::Matrix2d DDs[6];
    DDs[0] << 1, 0, 0, 0;   // x1
    DDs[1] << 0, 0, 1, 0;   // y1
    DDs[2] << 0, 1, 0, 0;   // x2
    DDs[3] << 0, 0, 0, 1;   // y2
    DDs[4] << -1, -1, 0, 0;
    DDs[5] << 0, 0, -1, -1;

    Eigen::Matrix2d DF[6];  // derivative of F with respect to x1,y1,x2,...
    for(int i=0;i<6;i++) DF[i] = DDs[i]*Dm_inv;

    double J = F.determinant();

    double log_J = log(J);
    strain_energy_density = (mu/2.0)*((F*FT).trace()-2.0)-mu*log_J+(lambda/2.0)*log_J*log_J;

    // First Piola - Kirchhoff stress tensor
    Eigen::Matrix2d P = F*mu + FinvT*(lambda*log_J-mu);

    // forces on nodes 1 and 2
    Eigen::Matrix2d H = -W*P*Dm_inv.transpose();

    // energy gradient with respect to x1,y1,x2,y2,x3,y3
    DE(0) = H(0,0);
    DE(1) = H(1,0);
    DE(2) = H(0,1);
    DE(3) = H(1,1);
    DE(4) = -DE(0)-DE(2);
    DE(5) = -DE(1)-DE(3);

    for(int i=0;i<6;i++)
    {
        Eigen::Matrix2d dP = mu*DF[i] + (mu-lambda*log_J)*FinvT*DF[i].transpose()*FinvT + lambda*(Finv*DF[i]).trace()*FinvT;
        Eigen::Matrix2d dH = -W*dP*Dm_inv.transpose();
        HE(0,i) = dH(0,0);
        HE(1,i) = dH(1,0);
        HE(2,i) = dH(0,1);
        HE(3,i) = dH(1,1);
        HE(4,i) = -dH(0,0)-dH(0,1);
        HE(5,i) = -dH(1,0)-dH(1,1);
    }

    for(int i=0;i<3;i++)
    {
        int row = nds[i]->eqId;
        Eigen::Vector2d locDE = -DE.block(i*2,0,2,1);
        eq.AddToC(row, locDE);
        for(int j=0;j<3;j++)
        {
            int col = nds[j]->eqId;
            Eigen::Matrix2d locHE = -HE.block(i*2,j*2,2,2);
            eq.AddToQ(row, col, locHE);
        }
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




