#ifndef ELEMENT123_H
#define ELEMENT123_H



#include <Eigen/Core>
#include <Eigen/Geometry>

//#include "parameters_sim.h"
#include "node.h"

namespace icy { class Element; class Node; }

class icy::Element
{
public:
    void Reset(void);

    icy::Node* nds[3];          // initialized when the geometry is loaded or remeshed
//    icy::Edge edges[3];        // element's edges opposite to nd 0,1,2
    icy::Element* adj_elems[3]; // nullptr if no adjacent element

    // at initial state
    double area_initial;
    void PrecomputeInitialArea();

    /*
    void UpdateSparseSystemEntries(LinearSystem &ls);
    void ComputeElasticForce(LinearSystem &ls, SimParams &prms, double timeStep,
                             Eigen::Matrix3d &elasticityMatrix,
                             Eigen::Matrix2d &D_mats);

    void ComputeMatrices(SimParams &prms,
                         Eigen::Matrix3d &elasticityMatrix,
                         Eigen::Matrix2d &D_mats,
                         Eigen::Matrix<double,3,DOFS*3> &bmat_b,
                         Eigen::Matrix<double,2,DOFS*3> (&bmat_s)[3],
                         Eigen::Matrix<double,3,DOFS*3> &bmat_m,
                         Eigen::Matrix<double,DOFS*3,DOFS*3> &K);
*/

};

#endif // ELEMENT123_H
