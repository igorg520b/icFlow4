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

}


void icy::Element::ComputeEquationEntries(EquationOfMotionSolver &eq, SimParams &prms, double timeStep)
{

}


/*
void icy::Element::UpdateSparseSystemEntries(LinearSystem &ls)
{
    if(nds[0]->lsId < 0 && nds[1]->lsId < 0 && nds[2]->lsId < 0) return;

    // reserve non-zero entries in the sparse matrix structure
    ls.AddElementToStructure(nds[0]->lsId, nds[1]->lsId);
    ls.AddElementToStructure(nds[0]->lsId, nds[2]->lsId);
    ls.AddElementToStructure(nds[1]->lsId, nds[2]->lsId);
}



void icy::Element::ComputeMatrices(SimParams &prms,
                     Eigen::Matrix3d &elasticityMatrix,
                     Eigen::Matrix2d &D_mats,
                     Eigen::Matrix<double,3,DOFS*3> &bmat_b,
                     Eigen::Matrix<double,2,DOFS*3> (&bmat_s)[3],
                     Eigen::Matrix<double,3,DOFS*3> &bmat_m,
                     Eigen::Matrix<double,DOFS*3,DOFS*3> &K)
{
    // translate the element
    Eigen::Vector3d p1, p2;
    p1 = nds[1]->x_initial - nds[0]->x_initial;
    p2 = nds[2]->x_initial - nds[0]->x_initial;

    normal_initial = p1.cross(p2);
    area_initial=normal_initial.norm()/2;
    if(area_initial<1e-8) {
        qDebug() << "element " << nds[0]->locId << ", " << nds[1]->locId << ", " << nds[2]->locId;
        qDebug() << "area" << area_initial;
        qDebug() << "nd0 " << nds[0]->x_initial.x() << ", " << nds[0]->x_initial.y();
        qDebug() << "nd1 " << nds[1]->x_initial.x() << ", " << nds[1]->x_initial.y();
        qDebug() << "nd2 " << nds[2]->x_initial.x() << ", " << nds[2]->x_initial.y();
        throw std::runtime_error("degenerate element created");
    }
    normal_initial.normalize();

    // (use 1-based indices, yij = yi-yj)
    double x1 = 0;
    double x2 = p1.x();
    double x3 = p2.x();
    double y1 = 0;
    double y2 = p1.y();
    double y3 = p2.y();

    double A2 = 2*area_initial;
    double y23, y31, y12, x32, x13, x21;
    y23 = y2-y3;
    y31 = y3-y1;
    y12 = y1-y2;
    x32 = x3-x2;
    x13 = x1-x3;
    x21 = x2-x1;

    // derivatives of shape functions
    double dN1dx = y23/A2;
    double dN2dx = y31/A2;
    double dN3dx = y12/A2;
    double dN1dy = x32/A2;
    double dN2dy = x13/A2;
    double dN3dy = x21/A2;

    // bending strain-displacement matrix
    bmat_b <<
         0,0,0,dN1dx,0,       0,0,0,dN2dx,0,        0,0,0,dN3dx,0,
         0,0,0,0,dN1dy,       0,0,0,0,dN2dy,        0,0,0,0,dN3dy,
         0,0,0,dN1dy,dN1dx,  0,0,0,dN2dy,dN2dx,   0,0,0,dN3dy,dN3dx;

    // membrane strain-displacement matrix
    bmat_m <<
         dN1dx,0,0,0,0,       dN2dx,0,0,0,0,      dN3dx,0,0,0,0,
         0,dN1dy,0,0,0,       0,dN2dy,0,0,0,      0,dN3dy,0,0,0,
         dN1dy,dN1dx,0,0,0,   dN2dy,dN2dx,0,0,0,  dN3dy,dN3dx,0,0,0;

    // shear
    for(int i=0;i<3;i++)
        bmat_s[i] <<
           0,0,dN1dx,-N[i][0],0,   0,0,dN2dx,-N[i][1],0,     0,0,dN3dx,-N[i][2], 0,
           0,0,dN1dy,0,-N[i][0],   0,0,dN2dy,0,-N[i][1],     0,0,dN3dy,0,-N[i][2];

    double thickness = prms.Thickness;
    // K and M depend on rho, Young's modulus and Poisson's ratio,
    // therefore they are computed after these parameters are set

//    Eigen::Matrix<double,DOFS*3,DOFS*3> K_b, K_m, K_s;
//    K = Eigen::Matrix<double,DOFS*3,DOFS*3>::Zero();
    K = bmat_m.transpose()*elasticityMatrix*bmat_m*(area_initial*thickness);

    double coeff = area_initial*thickness*thickness*thickness/12.0;
    K += bmat_b.transpose()*elasticityMatrix*bmat_b*coeff;

    for(int i=0;i<3;i++) K += bmat_s[i].transpose()*D_mats*bmat_s[i]*(thickness*area_initial/3.0);
}

void icy::Element::ComputeElasticForce(icy::LinearSystem &ls, icy::SimParams &prms, double,
                                       Eigen::Matrix3d &elasticityMatrix,
                                       Eigen::Matrix2d &D_mats)
{
    if(nds[0]->lsId < 0 && nds[1]->lsId < 0 && nds[2]->lsId < 0) return;

    Eigen::Matrix<double,3,DOFS*3> bmat_b;
    Eigen::Matrix<double,2,DOFS*3> bmat_s[3];   // 3 gauss points
    Eigen::Matrix<double,3,DOFS*3> bmat_m;
    Eigen::Matrix<double,DOFS*3,DOFS*3> K;    // element stiffness matrix (3 gauss points)
    ComputeMatrices(prms, elasticityMatrix, D_mats, bmat_b, bmat_s, bmat_m, K);

    Eigen::Matrix<double,DOFS*3,1> un;
    Eigen::Matrix<double,DOFS*3,1> ut;
    un << nds[0]->un, nds[1]->un, nds[2]->un;
    ut << nds[0]->ut, nds[1]->ut, nds[2]->ut;

    // calculate elastic forces and Hessian at step n+1

    // absolute position is not important, the origin is set at node 0
    Eigen::Matrix<double,DOFS*3,1> Fn, Fnp1;     // internal force at steps n and n+1
    Eigen::Matrix<double,DOFS*3,DOFS*3> &dFnp1=K;

//    FdF(un, Fn, nullptr);
//    FdF(ut, Fnp1, &dFnp1);
    Fn = K*un;
    Fnp1 = K*ut;

    // combine the results into a linearized equation of motion with HHT-alpha integration scheme
    double alpha = prms.HHTalpha;
    Eigen::Matrix<double,DOFS*3,1> F;      // right-hand side of the equation is equal to -F
    Eigen::Matrix<double,DOFS*3,DOFS*3> dF;

    F = Fn*alpha + Fnp1*(1-alpha);
    dF= dFnp1*(1-alpha);

#ifdef QT_DEBUG
    // assert
    for(int i=0;i<DOFS*3;i++)
        for(int j=0;j<DOFS*3;j++)
            if(std::isnan(dF(i,j)))
                throw std::runtime_error("elem.ComputeElasticForce: dF contains NaN");
#endif

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

