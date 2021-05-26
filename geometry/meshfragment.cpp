#include "meshfragment.h"
#include <gmsh.h>
#include <unordered_set>


void icy::MeshFragment::GenerateBrick(double ElementSize)
{
    deformable = true;

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

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", ElementSize);

    GetFromGmsh();
}



void icy::MeshFragment::GenerateIndenter(double ElementSize)
{
    deformable = false;

    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    double height = 1;
    double radius = 0.15;
//    int point1 = gmsh::model::occ::addPoint(0, height+radius*1.1, 0, 1.0);

    int ellipseTag = gmsh::model::occ::addEllipse(0, height+radius*1.1, 0, radius, radius/2);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.CharacteristicLengthMax", ElementSize);

    GetFromGmsh();
}



void icy::MeshFragment::GetFromGmsh()
{
    gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary",0);
    gmsh::model::mesh::generate(deformable ? 2 : 1);


    // retrieve the result
    elems.clear();
    nodes.clear();
    boundary_nodes.clear();
    boundary_edges.clear();
    freeNodeCount = 0;


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

    std::unordered_set<std::size_t> tagSet; // only keep nodes from the tris and edges
    for(std::size_t &tag : nodeTagsInTris) tagSet.insert(tag);
    for(std::size_t &tag : nodeTagsInEdges) tagSet.insert(tag);

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
        nd.locId = count;
        nd.x_initial << x, y;
        nd.intended_position = nd.xt = nd.xn = nd.x_initial;
        if(y==0 || !deformable) nd.pinned=true;
        else freeNodeCount++;
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
    boundary_edges.resize(nodeTagsInEdges.size()/2);
    for(std::size_t i=0;i<nodeTagsInEdges.size()/2;i++)
    {
        int idx1 = mtags[nodeTagsInEdges[i*2+0]];
        int idx2 = mtags[nodeTagsInEdges[i*2+1]];
        Node *nd1 = &nodes[idx1];
        Node *nd2 = &nodes[idx2];
        boundary_edges[i]=std::make_pair(nd1,nd2);
        set_nds.insert(idx1);
        set_nds.insert(idx2);
    }

    boundary_nodes.resize(set_nds.size());
    std::copy(set_nds.begin(), set_nds.end(), boundary_nodes.begin());
    gmsh::clear();
}

/*

void MeshFragment::GenerateSelfCollisionTest(double ElementSize)
{
    deformable = true;

    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    GetFromGmsh();
}

void MeshFragment::GenerateCircle(double x, double y, double r, double ElementSize)
{
    deformable = true;


    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");

    GetFromGmsh();

}


void MeshFragment::GenerateCup(double ElementSize)
{
    deformable = false;

    // invoke Gmsh
    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 1);
    gmsh::model::add("block1");


    GetFromGmsh();
}
*/