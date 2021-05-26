#ifndef MESHFRAGMENT_H
#define MESHFRAGMENT_H

#include "element.h"
#include <vector>

namespace icy { class MeshFragment; }

class icy::MeshFragment
{
public:

    bool deformable;
    std::vector<icy::Node> nodes;
    std::vector<icy::Element> elems;
    std::vector<int> boundary_nodes;
    std::vector<std::pair<icy::Node*,icy::Node*>> boundary_edges;
    unsigned freeNodeCount;

    void GenerateBrick(double ElementSize);
    void GenerateIndenter(double ElementSize);
//    void GenerateCup(double ElementSize);
//    void GenerateSelfCollisionTest(double ElementSize);
//    void GenerateCircle(double x, double y, double r, double ElementSize);

private:
    void GetFromGmsh();

};

#endif // MESHFRAGMENT_H
