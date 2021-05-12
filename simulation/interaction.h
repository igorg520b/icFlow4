#ifndef INTERACTION_H
#define INTERACTION_H

#include<Eigen/Core>
#include "equationofmotionsolver.h"
#include "parameters_sim.h"

namespace icy { class Interaction; class Node;}

class icy::Interaction
{
public:
    unsigned ndA_idx, ndB_idx, ndP_idx;
    Eigen::Vector2d A, B, P, D;
    double t, dist;
    Node *ndA, *ndB, *ndP;

    void AddToSparsityStructure(EquationOfMotionSolver &eq);
    void Evaluate(EquationOfMotionSolver &eq, SimParams &prms);

private:
    void static distance(Eigen::Vector2d (&p)[4], double &d, double &t, Eigen::Matrix<double,6,1> &Dd, Eigen::Matrix<double,6,6> &DDd);

    void static potential(double dHat, double d, Eigen::Matrix<double,6,1> &Dd, Eigen::Matrix<double,6,6> &DDd,
                          double &p, Eigen::Matrix<double,6,1> &Dp, Eigen::Matrix<double,6,6> &DDp);
};

#endif // INTERACTION_H
