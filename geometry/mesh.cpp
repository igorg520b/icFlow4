#include "mesh.h"
#include "model.h"
#include <numeric>
#include <algorithm>
#include <iterator>
#include <gmsh.h>


//#include <bits/stdc++.h>


void icy::Mesh::Reset(double CharacteristicLengthMax)
{
    elems.clear();
    nodes.clear();
    boundary.clear();

    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    double width = 2;
    double height = 1;
    int point1 = gmsh::model::occ::addPoint(-width/2, 0, 0, 1.0);
    int point2 = gmsh::model::occ::addPoint(-width/2, height, 0, 1.0);
    int point3 = gmsh::model::occ::addPoint(width/2, height, 0, 1.0);
    int point4 = gmsh::model::occ::addPoint(width/2, 0, 0, 1.0);

    int line1 = gmsh::model::occ::addLine(point1, point2);
    int line2 = gmsh::model::occ::addLine(point2, point3);
    int line3 = gmsh::model::occ::addLine(point3, point4);
    int line4 = gmsh::model::occ::addLine(point4, point1);

    std::vector<int> curveTags;
    curveTags.push_back(line1);
    curveTags.push_back(line2);
    curveTags.push_back(line3);
    curveTags.push_back(line4);
    int loopTag = gmsh::model::occ::addCurveLoop(curveTags);

    std::vector<int> loops;
    loops.push_back(loopTag);
    gmsh::model::occ::addPlaneSurface(loops);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", CharacteristicLengthMax);
    gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary",0);
    gmsh::model::mesh::generate(2);


    // retrieve the result

    // nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords);

    // elems
    std::vector<std::size_t> trisTags, nodeTagsInTris;
    gmsh::model::mesh::getElementsByType(2, trisTags, nodeTagsInTris);

    // boundary
    std::vector<std::size_t> edgeTags, nodeTagsInEdges;
    gmsh::model::mesh::getElementsByType(1, edgeTags, nodeTagsInEdges);

    // nodeTags, nodeCoords, nodeTagsInTris => nodes, elems
    std::map<std::size_t, int> nodeTagsMap1; // nodeTag -> its sequential position in nodeTag
    for(std::size_t i=0;i<nodeTags.size();i++) nodeTagsMap1[nodeTags[i]] = i;

    std::unordered_set<std::size_t> tagSet; // only keep nodes from the tris
    for(std::size_t &tag : nodeTagsInTris) tagSet.insert(tag);

    // set the size of the resulting nodes array
    nodes.resize(tagSet.size());

    std::map<std::size_t, int> mtags; // nodeTag -> sequential position
    int count = 0;

    for(const std::size_t &tag : tagSet)
    {
        int idx1 = nodeTagsMap1[tag];
        double x = nodeCoords[idx1*3+0];
        double y = nodeCoords[idx1*3+1];

        icy::Node &nd = nodes[count];
        nd.Reset();
        nd.id = count;
        nd.x_initial << x, y;
        nd.xt = nd.xn = nd.x_initial;
        if(y==0) nd.pinned=true;
        mtags[tag] = count;
        count++;
    }

    // resulting elements array
    elems.resize(nodeTagsInTris.size()/3);

    for(std::size_t i=0;i<nodeTagsInTris.size()/3;i++)
    {
        icy::Element &elem = elems[i];
        for(int j=0;j<3;j++) elem.nds[j] = &(nodes[mtags[nodeTagsInTris[i*3+j]]]);
        elem.PrecomputeInitialArea();
        if(elem.area_initial == 0) throw std::runtime_error("icy::Mesh::Reset - element's area is zero");
        if(elem.area_initial < 0)
        {
            for(int j=0;j<3;j++) elem.nds[2-j] = &(nodes[mtags[nodeTagsInTris[i*3+j]]]);
            elem.PrecomputeInitialArea();
            if(elem.area_initial < 0) throw std::runtime_error("icy::Mesh::Reset - error");
        }
        for(int j=0;j<3;j++) elem.nds[j]->area += elem.area_initial/3;
    }

    std::unordered_set<int> set_nds; // set of boundary nodes
    boundary.resize(nodeTagsInEdges.size()/2);
    for(std::size_t i=0;i<nodeTagsInEdges.size()/2;i++)
    {
        int idx1 = mtags[nodeTagsInEdges[i*2+0]];
        int idx2 = mtags[nodeTagsInEdges[i*2+1]];
        boundary[i]=std::make_pair(idx1,idx2);
        set_nds.insert(idx1);
        set_nds.insert(idx2);
    }

    deformable_boundary_nodes.resize(set_nds.size());
    std::copy(set_nds.begin(), set_nds.end(), deformable_boundary_nodes.begin());

    GenrateIndenter(CharacteristicLengthMax);

}


void icy::Mesh::GenrateIndenter(double CharacteristicLengthMax)
{
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    double height = 1;
    double radius = 0.15;
//    int point1 = gmsh::model::occ::addPoint(0, height+radius*1.1, 0, 1.0);

    int ellipseTag = gmsh::model::occ::addEllipse(0, height+radius*1.1, 0, radius, radius/2);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", CharacteristicLengthMax);
    gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary",0);
    gmsh::model::mesh::generate(2);

    // nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords);
//    std::cout << "nds " << nodeTags.size() << std::endl;

    std::map<std::size_t, int> nodeTagsMap1; // nodeTag -> its sequential position in nodeTag
    for(std::size_t i=0;i<nodeTags.size();i++) nodeTagsMap1[nodeTags[i]] = i;

    // boundary
    std::vector<std::size_t> edgeTags, nodeTagsInEdges;
    gmsh::model::mesh::getElementsByType(1, edgeTags, nodeTagsInEdges);
    std::unordered_set<unsigned> boundary_nds_idxs;
    for(unsigned i=0;i<edgeTags.size();i++)
    {
        boundary_nds_idxs.insert(nodeTagsInEdges[i*2+0]);
        boundary_nds_idxs.insert(nodeTagsInEdges[i*2+1]);
    }


    std::map<std::size_t, int> nodeTagsMap2; // nodeTag -> its sequential position in nodes_indenter
    unsigned count = 0;
    nodes_indenter.resize(boundary_nds_idxs.size());
    for(const unsigned tag : boundary_nds_idxs)
    {
        int idx1 = nodeTagsMap1[tag];
        double x = nodeCoords[idx1*3+0];
        double y = nodeCoords[idx1*3+1];

        nodeTagsMap2[tag]=count;

        icy::Node &nd = nodes_indenter[count];
        nd.Reset();
        nd.id = count;
        nd.x_initial << x, y;
        nd.xt = nd.xn = nd.x_initial;
        nd.pinned=true;

        count++;
    }

    boundary_indenter.clear();
    for(unsigned i=0;i<edgeTags.size();i++)
    {
        int idx1 = nodeTagsInEdges[i*2+0];
        int idx2 = nodeTagsInEdges[i*2+1];
        boundary_indenter.push_back(std::make_pair(nodeTagsMap2[idx1],nodeTagsMap2[idx2]));
    }

}


void icy::Mesh::DetectContactPairs(double distance_threshold)
{

    indenter_boundary_vs_deformable_nodes.clear();
    deformable_boundary_vs_indenter_nodes.clear();

    for(unsigned b_idx=0; b_idx<boundary_indenter.size(); b_idx++)
        for(unsigned n_idx=0; n_idx<deformable_boundary_nodes.size(); n_idx++)
        {
            Interaction i;
            i.ndA_idx = boundary_indenter[b_idx].first;
            i.ndB_idx = boundary_indenter[b_idx].second;
            i.ndP_idx = deformable_boundary_nodes[n_idx];

            icy::Node &ndA = nodes_indenter[i.ndA_idx];
            icy::Node &ndB = nodes_indenter[i.ndB_idx];
            icy::Node &ndP = nodes[i.ndP_idx];

            i.dist = SegmentPointDistance(ndA.xt, ndB.xt, ndP.xt, i.t);
            if(i.dist <= distance_threshold)
                indenter_boundary_vs_deformable_nodes.push_back(i);
        }

    qDebug() << "icy::Mesh::DetectContactPairs(): ib_dn " << indenter_boundary_vs_deformable_nodes.size();
}

double icy::Mesh::SegmentPointDistance(Eigen::Vector2d A, Eigen::Vector2d B, Eigen::Vector2d P, double &t)
{
    Eigen::Vector2d seg = B-A;
    Eigen::Vector2d v = P-A;
    t = v.dot(seg)/seg.squaredNorm();
    t = std::clamp(t, 0.0, 1.0);
    Eigen::Vector2d D = A+seg*t;
    double dist = (D-P).norm();
    return dist;
}
