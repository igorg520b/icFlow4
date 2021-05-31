#include "mesh.h"
#include "model.h"
#include <numeric>
#include <algorithm>
#include <iterator>





icy::Mesh::Mesh()
{
    InitializeLUT(2);
    visualized_values->SetName("visualized_values");

    ugrid_deformable->SetPoints(points_deformable);
    dataSetMapper_deformable->SetInputData(ugrid_deformable);
    dataSetMapper_deformable->UseLookupTableScalarRangeOn();
    dataSetMapper_deformable->SetLookupTable(hueLut);

    actor_mesh_deformable->SetMapper(dataSetMapper_deformable);
    actor_mesh_deformable->GetProperty()->VertexVisibilityOff();
    actor_mesh_deformable->GetProperty()->EdgeVisibilityOn();
    actor_mesh_deformable->GetProperty()->SetColor(0.84938, 0.872213, 0.848103);
    actor_mesh_deformable->GetProperty()->SetEdgeColor(90.0/255.0, 90.0/255.0, 97.0/255.0);
    actor_mesh_deformable->GetProperty()->LightingOff();
    actor_mesh_deformable->GetProperty()->ShadingOff();
    actor_mesh_deformable->GetProperty()->SetInterpolationToFlat();
    actor_mesh_deformable->PickableOff();
    actor_mesh_deformable->GetProperty()->SetLineWidth(1);

    // boundary that encompasses all objects
    ugrid_boundary_all->SetPoints(points_deformable);
    dataSetMapper_boundary_all->SetInputData(ugrid_boundary_all);
    actor_boundary_all->SetMapper(dataSetMapper_boundary_all);
    actor_boundary_all->GetProperty()->EdgeVisibilityOn();
    actor_boundary_all->GetProperty()->VertexVisibilityOn();
    actor_boundary_all->GetProperty()->SetColor(0,0,0.4);
    actor_boundary_all->GetProperty()->SetEdgeColor(0,0,0.4);
    actor_boundary_all->GetProperty()->SetVertexColor(0.4,0,0);
    actor_boundary_all->GetProperty()->SetPointSize(5);
    actor_boundary_all->GetProperty()->SetLineWidth(3);
    actor_boundary_all->PickableOff();

    // indenter's intended position
    ugrid_indenter_intended->SetPoints(points_indenter_intended);
    dataSetMapper_indenter_intended->SetInputData(ugrid_indenter_intended);
    actor_boundary_intended_indenter->SetMapper(dataSetMapper_indenter_intended);
    actor_boundary_intended_indenter->GetProperty()->EdgeVisibilityOn();
    actor_boundary_intended_indenter->GetProperty()->VertexVisibilityOn();
    actor_boundary_intended_indenter->GetProperty()->SetColor(0.3,0,0.4);
    actor_boundary_intended_indenter->GetProperty()->SetEdgeColor(0.6,0,0.04);
    actor_boundary_intended_indenter->GetProperty()->SetVertexColor(0.04,0,0);
    actor_boundary_intended_indenter->GetProperty()->SetPointSize(2);
    actor_boundary_intended_indenter->GetProperty()->SetLineWidth(1);

    // collisions
    ugrid_collisions->SetPoints(points_collisions);
    mapper_collisions->SetInputData(ugrid_collisions);
    actor_collisions->SetMapper(mapper_collisions);
    actor_collisions->GetProperty()->EdgeVisibilityOn();
    actor_collisions->GetProperty()->VertexVisibilityOn();
    actor_collisions->GetProperty()->SetColor(0.3,0,0.4);
    actor_collisions->GetProperty()->SetEdgeColor(0.6,0,0.04);
    actor_collisions->GetProperty()->SetVertexColor(0.04,0,0);
    actor_collisions->GetProperty()->SetPointSize(2);
    actor_collisions->GetProperty()->SetLineWidth(1);

    collision_interactions.reserve(10000);
    broadphase_list.reserve(100000);

    root_contact.isLeaf = root_ccd.isLeaf = false;
    root_contact.test_self_collision = root_ccd.test_self_collision = true;
}


void icy::Mesh::Reset(double CharacteristicLengthMax)
{
    brick.GenerateBrick(CharacteristicLengthMax);
    indenter.GenerateIndenter(CharacteristicLengthMax);
    allMeshes.clear();
    allMeshes.push_back(&brick);
    allMeshes.push_back(&indenter);
    RegenerateVisualizedGeometry();
    tree_update_counter=0;
}

void icy::Mesh::RegenerateVisualizedGeometry()
{
    allNodes.clear();
    allElems.clear();
    allBoundaryEdges.clear();
    unsigned count = 0;
    freeNodeCount = 0;

    global_leafs_ccd.clear();
    global_leafs_contact.clear();
    fragmentRoots_ccd.clear();
    fragmentRoots_contact.clear();

    for(MeshFragment *mf : allMeshes)
    {
        for(unsigned i=0;i<mf->nodes.size();i++)
        {
            Node *nd = &mf->nodes[i];
            nd->globId = count++;
            if(nd->pinned) nd->eqId=-1;
            else nd->eqId=freeNodeCount++;
            allNodes.push_back(nd);
        }
        mf->GenerateLeafs();

        for(unsigned i=0;i<mf->elems.size();i++) allElems.push_back(&mf->elems[i]);
        allBoundaryEdges.insert(allBoundaryEdges.end(), mf->boundary_edges.begin(), mf->boundary_edges.end());
        global_leafs_ccd.insert(global_leafs_ccd.end(), mf->leafs_for_ccd.begin(), mf->leafs_for_ccd.end());
        global_leafs_contact.insert(global_leafs_contact.end(),mf->leafs_for_contact.begin(), mf->leafs_for_contact.end());
        fragmentRoots_ccd.push_back(&mf->root_ccd);
        fragmentRoots_contact.push_back(&mf->root_contact);
    }

    points_deformable->SetNumberOfPoints(allNodes.size());
    cellArray_deformable->Reset();

    // deformable material - elements
    for(icy::Element *tr : allElems)
    {
        vtkIdType pts[3] = {tr->nds[0]->globId, tr->nds[1]->globId, tr->nds[2]->globId};
        cellArray_deformable->InsertNextCell(3, pts);
    }
    ugrid_deformable->SetCells(VTK_TRIANGLE, cellArray_deformable);


    // all boundaries
    cellArray_boundary_all->Reset();
    for(auto edge : allBoundaryEdges)
    {
        vtkIdType pts[2] = {edge.first->globId, edge.second->globId};
        cellArray_boundary_all->InsertNextCell(2, pts);
    }
    ugrid_boundary_all->SetCells(VTK_LINE, cellArray_boundary_all);


    // intended position of the indenter
    points_indenter_intended->SetNumberOfPoints(indenter.nodes.size());
    cellArray_indenter_intended->Reset();
    for(auto edge : indenter.boundary_edges)
    {
        vtkIdType pts[2] = {edge.first->locId, edge.second->locId};
        cellArray_indenter_intended->InsertNextCell(2, pts);
    }
    ugrid_indenter_intended->SetCells(VTK_LINE, cellArray_indenter_intended);

}


void icy::Mesh::UpdateTree(float distance_threshold)
{
    qDebug() << "icy::Mesh::UpdateTree " << distance_threshold;
    // update leafs
    unsigned nLeafs = global_leafs_ccd.size();
#pragma omp parallel for
    for(unsigned i=0;i<nLeafs;i++)
    {
        BVHN *leaf_ccd = global_leafs_ccd[i];
        auto [nd1Idx,nd2Idx] = leaf_ccd->feature;
        Node *nd1 = allNodes[nd1Idx];
        Node *nd2 = allNodes[nd2Idx];
        kDOP8 &box_ccd = leaf_ccd->box;
        box_ccd.Reset();
        box_ccd.Expand(nd1->xn.x(), nd1->xn.y());
        box_ccd.Expand(nd2->xn.x(), nd2->xn.y());
        box_ccd.Expand(nd1->xt.x(), nd1->xt.y());
        box_ccd.Expand(nd2->xt.x(), nd2->xt.y());

        BVHN *leaf_contact = global_leafs_contact[i];
        nd1 = allNodes[leaf_contact->feature.first];
        nd2 = allNodes[leaf_contact->feature.second];

        kDOP8 &box_contact = leaf_contact->box;
        box_contact.Reset();
        box_contact.Expand(nd1->xt.x(), nd1->xt.y());
        box_contact.Expand(nd2->xt.x(), nd2->xt.y());
        box_contact.ExpandBy(distance_threshold);
    }

    qDebug() << "leafs updated " << nLeafs;

    // update or build the rest of the tree
    // TODO: parallel
    if(tree_update_counter%10 != 0)
    {
        qDebug() << "updating tree " << tree_update_counter%10;
        root_ccd.Update();
        root_contact.Update();
    }
    else
    {
        qDebug() << "building tree";
        BVHN::BVHNFactory.releaseAll(); // does not release root nodes of mesh fragments
        for(MeshFragment *mf : allMeshes)
        {
            mf->root_ccd.Build(&mf->leafs_for_ccd,0);
            mf->root_contact.Build(&mf->leafs_for_contact,0);
        }
        root_ccd.Build(&fragmentRoots_ccd,0);
        root_contact.Build(&fragmentRoots_contact,0);
    }
    tree_update_counter++;
    qDebug() << "icy::Mesh::UpdateTree " << tree_update_counter;
}






void icy::Mesh::ChangeVisualizationOption(int option)
{
    qDebug() << "icy::Model::ChangeVisualizationOption " << option;
    if(VisualizingVariable == option) return; // option did not change
    VisualizingVariable = option;

    if(VisualizingVariable == (int)icy::Model::VisOpt::none)
    {
        dataSetMapper_deformable->ScalarVisibilityOff();
        ugrid_deformable->GetPointData()->RemoveArray("visualized_values");
        ugrid_deformable->GetCellData()->RemoveArray("visualized_values");
        return;
    }
    else
    {
        ugrid_deformable->GetPointData()->RemoveArray("visualized_values");
        ugrid_deformable->GetCellData()->AddArray(visualized_values);
        ugrid_deformable->GetCellData()->SetActiveScalars("visualized_values");
        dataSetMapper_deformable->SetScalarModeToUseCellData();
        dataSetMapper_deformable->ScalarVisibilityOn();
    }
    UpdateValues();
}






void icy::Mesh::UnsafeUpdateGeometry()
{
    for(icy::Node *nd : allNodes) points_deformable->SetPoint((vtkIdType)nd->globId, nd->xn.data());
    points_deformable->Modified();

    // indenter intended points
    for(icy::Node &nd : indenter.nodes)
        points_indenter_intended->SetPoint((vtkIdType)nd.locId, nd.intended_position.data());
    points_indenter_intended->Modified();

    if(VisualizingVariable != icy::Model::VisOpt::none) UpdateValues();


    // collisions
    points_collisions->SetNumberOfPoints(collision_interactions.size()*2);
    cellArray_collisions->Reset();
    for(unsigned i=0;i<collision_interactions.size();i++)
    {
        vtkIdType pts[2] = {2*i, 2*i+1};
        cellArray_collisions->InsertNextCell(2, pts);
        Interaction &intr = collision_interactions[i];
        points_collisions->SetPoint((vtkIdType)2*i, intr.ndP->xt.data());
        points_collisions->SetPoint((vtkIdType)2*i+1, intr.D.data());
    }
    points_indenter_intended->Modified();
    ugrid_collisions->SetCells(VTK_LINE, cellArray_collisions);
    actor_collisions->Modified();



/*
    glyph_int_data->SetNumberOfValues(mesh.nodes.size());

    for(std::size_t i=0;i<mesh.nodes.size();i++)
    {
        int value = 0;
        if(mesh.nodes[i].pinned)
            value = mesh.nodes[i].selected ? 2 : 1;
        glyph_int_data->SetValue(i,value);
    }
    glyph_hueLut->SetTableRange(-0.5,5.5);

    glyph_mapper->SetScalarModeToUsePointData();
    poly_data->GetPointData()->AddArray(glyph_int_data);
    poly_data->GetPointData()->SetActiveScalars("glyph_int_data");
*/

}

void icy::Mesh::AddToNarrowListIfNeeded(Node *ndA, Node *ndB, Node *ndP, double distance_threshold)
{
    Eigen::Vector2d D;
    double dist = icy::Interaction::SegmentPointDistance(ndA->xt, ndB->xt, ndP->xt, D);
    if(dist < distance_threshold)
    {
        Interaction i;
        i.ndA = ndA;
        i.ndB = ndB;
        i.ndP = ndP;
        i.D = D;
        collision_interactions.push_back(i);
    }
}


void icy::Mesh::DetectContactPairs(double distance_threshold)
{
    unsigned nEdges = allBoundaryEdges.size();
    broadphase_list.clear();
    for(unsigned idx1=0;idx1<nEdges;idx1++)
        for(unsigned idx2=0;idx2<idx1;idx2++)
            broadphase_list.push_back(std::make_pair(idx1,idx2));

    unsigned nBroadList = broadphase_list.size();

    collision_interactions.clear();
#pragma omp parallel for
    for(unsigned i=0;i<nBroadList;i++)
    {
        auto [idx1,idx2] = broadphase_list[i];
        auto [nd1, nd2] = allBoundaryEdges[idx1];
        auto [nd3, nd4] = allBoundaryEdges[idx2];

        AddToNarrowListIfNeeded(nd1, nd2, nd3, distance_threshold);
        AddToNarrowListIfNeeded(nd1, nd2, nd4, distance_threshold);
        AddToNarrowListIfNeeded(nd3, nd4, nd1, distance_threshold);
        AddToNarrowListIfNeeded(nd3, nd4, nd2, distance_threshold);
    }
}

void icy::Mesh::UpdateValues()
{
    if(allNodes.size()==0)
    {
        dataSetMapper_deformable->ScalarVisibilityOff();
        ugrid_deformable->GetPointData()->RemoveArray("visualized_values");
        ugrid_deformable->GetCellData()->RemoveArray("visualized_values");
        return;
    }

    switch((icy::Model::VisOpt)VisualizingVariable)
    {
        case icy::Model::VisOpt::elem_area:
        visualized_values->SetNumberOfValues(allElems.size());
        for(size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->area_initial);
        break;

    case icy::Model::VisOpt::energy_density:
        visualized_values->SetNumberOfValues(allElems.size());
        for(size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->strain_energy_density);
        break;

    case icy::Model::VisOpt::stress_xx:
        visualized_values->SetNumberOfValues(allElems.size());
        for(size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->CauchyStress(0,0));
        break;

    case icy::Model::VisOpt::stress_yy:
        visualized_values->SetNumberOfValues(allElems.size());
        for(size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->CauchyStress(1,1));
        break;

    case icy::Model::VisOpt::stress_hydrostatic:
        visualized_values->SetNumberOfValues(allElems.size());
        for(size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->CauchyStress.trace()/2);
        break;

    case icy::Model::VisOpt::non_symm_measure:
        visualized_values->SetNumberOfValues(allElems.size());
        for(size_t i=0;i<allElems.size();i++)
        {
            Element *elem = allElems[i];
            double value1 = elem->CauchyStress(0,1)-elem->CauchyStress(1,0);
            double value = value1*value1;
            visualized_values->SetValue(i, value);
        }
        break;

    case icy::Model::VisOpt::ps1:
        visualized_values->SetNumberOfValues(allElems.size());
        for(size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->principal_stress1);
        break;

    case icy::Model::VisOpt::ps2:
        visualized_values->SetNumberOfValues(allElems.size());
        for(size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->principal_stress2);
        break;

    case icy::Model::VisOpt::shear_stress:
        visualized_values->SetNumberOfValues(allElems.size());
        for(size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->max_shear_stress);
        break;

    case icy::Model::VisOpt::volume_change:
        visualized_values->SetNumberOfValues(allElems.size());
        for(size_t i=0;i<allElems.size();i++) visualized_values->SetValue(i, allElems[i]->volume_change);
        break;

    default:
        break;
    }

    visualized_values->Modified();

    double minmax[2];
    visualized_values->GetValueRange(minmax);
    hueLut->SetTableRange(minmax[0], minmax[1]);
    if(VisualizingVariable == icy::Model::VisOpt::volume_change)
    {
        double range = minmax[1]-minmax[0];
        hueLut->SetTableRange(1-range*0.75,1+range*0.75);
    }
}

void icy::Mesh::InitializeLUT(int table)
{
    const int n = 51;
    hueLut->SetNumberOfTableValues(n);

    if(table==1)
    for ( int i=0; i<n; i++)
            hueLut->SetTableValue(i, (double)lutArrayTerrain[i][0],
                    (double)lutArrayTerrain[i][1],
                    (double)lutArrayTerrain[i][2], 1.0);
    else if(table==2)
        for ( int i=0; i<n; i++)
                hueLut->SetTableValue(i, (double)lutArrayTemperatureAdj[i][0],
                        (double)lutArrayTemperatureAdj[i][1],
                        (double)lutArrayTemperatureAdj[i][2], 1.0);
}
